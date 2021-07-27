#include <ros/ros.h>
#include <ros/package.h>

#include <armadillo>

#include <gmp_lib/gmp_lib.h>

#include <gmp_test/utils/utils.h>

using namespace as64_;

// ==================================
// ------------   MAIN  -------------
// ==================================

arma::vec maxVec(const arma::vec &v1, const arma::vec &v2)
{
  arma::vec v = v1;
  for (int i=0; i<v.size(); i++)
  { if (v(i) < v2(i)) v(i) = v2(i); }
  return v;
}

arma::vec minVec(const arma::vec &v1, const arma::vec &v2)
{
  arma::vec v = v1;
  for (int i=0; i<v.size(); i++)
  { if (v(i) > v2(i)) v(i) = v2(i); }
  return v;
}

int main(int argc, char **argv)
{
  ros::init(argc, argv, "optimize_gmp_test_node");

  std::string in_filename = ros::package::getPath(PROJECT_NAME_) + "/src/optimization/constr_opt/data/pos_data.bin";
  std::string constr_filename = ros::package::getPath(PROJECT_NAME_) + "/src/optimization/constr_opt/data/constraints.bin";
  std::string out_filename = ros::package::getPath(PROJECT_NAME_) + "/src/optimization/constr_opt/data/results.bin";

  // ===========  Load training data  ===========
  arma::rowvec Timed;
  arma::mat Pd_data, dPd_data, ddPd_data;
  gmp_::FileIO fid(in_filename, gmp_::FileIO::in);
  fid.read("Timed", Timed);
  fid.read("Pd_data", Pd_data);
  fid.read("dPd_data", dPd_data);
  fid.read("ddPd_data", ddPd_data);
  fid.close();

  // ===========  initialize and train GMP  ===========
  std::string train_method = "LS";
  int N_kernels = 30;
  double kernels_std_scaling = 1.5;
  int n_dof = Pd_data.n_rows;
  gmp_::GMP::Ptr gmp( new gmp_::GMP(n_dof, N_kernels, kernels_std_scaling) );
  Timer::tic();
  arma::vec offline_train_mse;
  gmp->train(train_method, Timed/Timed.back(), Pd_data, &offline_train_mse);
  std::cout << "offline_train_mse = \n" << offline_train_mse << "\n";
  Timer::toc();

  // gmp->setScaleMethod(TrajScale.ROT_MIN_SCALE);

  double taud = Timed.back();
  arma::vec yd0 = gmp->getYd(0); // Pd_data(1);
  arma::vec ygd = gmp->getYd(1); // Pd_data(end);

  double kt = 1.0; // temporal scaling
  arma::mat ks = arma::diagmat( arma::vec({1.3, 1.4, 1.5}) ); // spatial scaling
  double tau = taud/kt;
  arma::vec y0 = yd0 + 0.15;
  arma::vec yg = ks*(ygd - yd0) + y0;
  gmp->setY0(y0);
  gmp->setGoal(yg);

  // ===========  calculate scaled demo trajectory  ===========
  arma::mat T_sc = gmp->getScaling();
  arma::rowvec Time2 = Timed / kt;
  arma::mat P_data2 = T_sc*Pd_data + arma::repmat(-T_sc*Pd_data.col(0) + y0, 1, Pd_data.n_cols);
  arma::mat dP_data2 = kt*T_sc*dPd_data;
  arma::mat ddP_data2 = std::pow(kt,2)*T_sc*ddPd_data;

  // ===========  Load constraints  ===========
  arma::rowvec x_pos_lim, xeq_pos;
  arma::mat pos_lb, pos_ub, pos_eq;

  arma::rowvec x_vel_lim, xeq_vel;
  arma::mat vel_lb, vel_ub, vel_eq;

  arma::rowvec x_accel_lim, xeq_accel;
  arma::mat accel_lb, accel_ub, accel_eq;
  {
    gmp_::FileIO fid(constr_filename, gmp_::FileIO::in);
    // -----------------------------
    fid.read("t_pos_lim", x_pos_lim);
    x_pos_lim /= tau;
    fid.read("pos_lb", pos_lb);
    fid.read("pos_ub", pos_ub);
    fid.read("teq_pos", xeq_pos);
    xeq_pos /= tau;
    fid.read("pos_eq", pos_eq);
    fid.read("t_vel_lim", x_vel_lim);
    x_vel_lim /= tau;
    fid.read("vel_lb", vel_lb);
    fid.read("vel_ub", vel_ub);
    fid.read("teq_vel", xeq_vel);
    xeq_vel /= tau;
    fid.read("vel_eq", vel_eq);
    fid.read("t_accel_lim", x_accel_lim);
    x_accel_lim /= tau;
    fid.read("accel_lb", accel_lb);
    fid.read("accel_ub", accel_ub);
    fid.read("teq_accel", xeq_accel);
    xeq_accel /= tau;
    fid.read("accel_eq", accel_eq);
    // -----------------------------
    fid.close();
  }

  // ===========  Solve constr-optimization problem  ===========
  Timer::tic();

  gmp_::GMP_Opt gmp_opt(gmp.get());
  gmp_opt.setProblemOptions(true, true, false, 0.1, 1, 0.1);
  gmp_opt.setMotionDuration(tau);
  gmp_opt.setPosConstr(x_pos_lim, pos_lb, pos_ub, xeq_pos, pos_eq);
  gmp_opt.setVelConstr(x_vel_lim, vel_lb, vel_ub, xeq_vel, vel_eq);
  gmp_opt.setAccelConstr(x_accel_lim, accel_lb, accel_ub, xeq_accel, accel_eq);
  gmp_opt.optimize(200);
  std::cout << gmp_opt.getExitMsg() << "\n";

  Timer::toc();

  // ===========  Generate trajectory of the constr-optimized GMP  ===========
  arma::rowvec Time;
  arma::mat P_data;
  arma::mat dP_data;
  arma::mat ddP_data;

  arma::vec p = y0;
  arma::vec p_dot;
  arma::vec p_ddot;

  double t = 0;
  double dt = 0.01;

  while (true)
  {
    double x = t/tau;
    double x_dot = 1/tau;
    p = gmp->getYd(x);
    p_dot = gmp->getYdDot(x, x_dot);
    p_ddot = gmp->getYdDDot(x, x_dot, 0);

    Time = arma::join_horiz(Time, arma::vec{t});
    P_data = arma::join_horiz(P_data, p);
    dP_data = arma::join_horiz(dP_data, p_dot);
    ddP_data = arma::join_horiz(ddP_data, p_ddot);

    t = t + dt;

    if (t >= 1.05*tau) break;
  }

  {
    gmp_::FileIO fid(out_filename, gmp_::FileIO::out | gmp_::FileIO::trunc);
    // -----------------------------
    fid.write("Time", Time);
    fid.write("P_data", P_data);
    fid.write("dP_data", dP_data);
    fid.write("ddP_data", ddP_data);
    // -----------------------------
    fid.write("Time2", Time2);
    fid.write("P_data2", P_data2);
    fid.write("dP_data2", dP_data2);
    fid.write("ddP_data2", ddP_data2);
    // -----------------------------
    fid.close();
  }


}
