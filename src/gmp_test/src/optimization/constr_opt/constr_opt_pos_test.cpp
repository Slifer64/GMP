#include <ros/ros.h>
#include <ros/package.h>

#include <armadillo>

#include <gmp_lib/gmp_lib.h>

#include <gmp_test/utils/utils.h>

using namespace as64_;

void getGMPTrajectory(gmp_::GMP::Ptr gmp, double tau, const arma::vec &y0, const arma::vec &yg,
                arma::rowvec &Time, arma::mat &P_data, arma::mat &dP_data, arma::mat &ddP_data);

void getOptGMPTrajectory(gmp_::GMP::Ptr gmp, double tau, const arma::vec &y0, const arma::vec &yg,
        const arma::mat &pos_lim, const arma::vec &vel_lim, const arma::vec &accel_lim, bool opt_pos, bool opt_vel,
        arma::rowvec &Time, arma::mat &P_data, arma::mat &dP_data, arma::mat &ddP_data);

// ==================================
// ------------   MAIN  -------------
// ==================================

class Data
{
public:
  Data(const arma::rowvec &Time, const arma::mat &Pos, const arma::mat &Vel, const arma::mat &Accel,
  const std::string &linestyle, const std::vector<double> &color, const std::string &legend, bool plot3D, bool plot2D)
  {
    this->Time = Time;
    this->Pos = Pos;
    this->Vel = Vel;
    this->Accel = Accel;
    this->linestyle = linestyle;
    this->color = color;
    this->legend = legend;
    this->plot3D = plot3D;
    this->plot2D = plot2D;
  }

  void write(gmp_::FileIO &fid, const std::string &prefix) const
  {
    fid.write(prefix + "Time", Time);
    fid.write(prefix + "Pos", Pos);
    fid.write(prefix + "Vel", Vel);
    fid.write(prefix + "Accel", Accel);
    fid.write(prefix + "linestyle", linestyle);
    fid.write(prefix + "color", color);
    fid.write(prefix + "legend", legend);
    fid.write(prefix + "plot3D", plot3D);
    fid.write(prefix + "plot2D", plot2D);
  }

  arma::rowvec Time;
  arma::mat Pos, Vel, Accel;
  std::string linestyle;
  std::vector<double> color;
  std::string legend;
  bool plot3D, plot2D;
};


int main(int argc, char **argv)
{
  ros::init(argc, argv, "optimize_gmp_test_node");

  int case_ = 1;
  if (argc > 1) case_ = atoi(argv[1]);

  std::string in_filename = ros::package::getPath(PROJECT_NAME_) + "/src/optimization/constr_opt/data/constr_opt_pos_test_train_data.bin";
  std::string out_filename = ros::package::getPath(PROJECT_NAME_) + "/src/optimization/constr_opt/data/results.bin";

  // ===========  Load training data  ===========
  arma::rowvec Timed;
  arma::mat Pd_data, dPd_data, ddPd_data;
  {
    gmp_::FileIO fid(in_filename, gmp_::FileIO::in);
    fid.read("Timed", Timed);
    fid.read("Pd_data", Pd_data);
    fid.read("dPd_data", dPd_data);
    fid.read("ddPd_data", ddPd_data);
    fid.close();
  }

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
  std::cout << "Training finished: "; Timer::toc();

  double taud = Timed.back();
  arma::vec yd0 = gmp->getYd(0); // Pd_data.col(0);
  arma::vec ygd = gmp->getYd(1); // Pd_data(end);

  double kt = 1.5; // temporal scaling
  double tau = taud/kt;
  arma::vec y0 = yd0;

  arma::vec yg = ygd;
  arma::rowvec view_;

  if (case_ == 1)
  {
    yg += arma::vec( {0.1, -0.1, 0.2} );
    view_ = {171.5301, -2.363};
  }
  else
  {
    yg += arma::vec( {0.7, -0.7, 0.05} );
    view_ = {171.9421, -3.069};
  }

  // ======== Limits ==========
  arma::mat pos_lim = arma::mat( { {-1.2, -1.2, 0.2}, {1.2, 1.2, 0.6} } ).t();
  arma::vec vel_lim = {-0.3, 0.3};
  arma::vec accel_lim = {-0.4, 0.4};

  // ======== Generate trajectories ==========

  arma::rowvec Time;
  arma::mat P_data, dP_data, ddP_data;

  gmp_::FileIO fid(out_filename, gmp_::FileIO::out | gmp_::FileIO::trunc);

  // --------- Proportional scaling -----------
  gmp->setScaleMethod( gmp_::TrajScale::Ptr( new gmp_::TrajScale_Prop(3) ) );
  getGMPTrajectory(gmp, tau, y0, yg, Time, P_data, dP_data, ddP_data);

  Data(Time, P_data, dP_data, ddP_data, ":", {0, 0, 1}, "prop", true, true).write(fid, "1_");

  // --------- Rotational scaling -----------
  gmp_::TrajScale::Ptr traj_sc( new gmp_::TrajScale_Rot_wb() );
  dynamic_cast<gmp_::TrajScale_Rot_wb *>(traj_sc.get())->setWorkBenchNormal(arma::vec({0,0,1}));
  gmp->setScaleMethod(traj_sc);

  getGMPTrajectory(gmp, tau, y0, yg, Time, P_data, dP_data, ddP_data);
  Data(Time, P_data, dP_data, ddP_data, ":", {0, 1, 1}, "rot-wb", true, true).write(fid, "2_");

  // --------- Demo -----------
  Data(Timed/kt, Pd_data, dPd_data*kt, ddPd_data*std::pow(kt,2), ":", {0.7, 0, 0}, "demo", true, false).write(fid, "3_");


  // --------- Optimized DMP -> VEL -----------
  getOptGMPTrajectory(gmp, tau, y0, yg, pos_lim, vel_lim, accel_lim, false, true, Time, P_data, dP_data, ddP_data);
  Data(Time, P_data, dP_data, ddP_data, "-", {0, 1, 0}, "opt-vel", true, true).write(fid, "4_");

  // // --------- Optimized DMP -> POS -----------
  // getOptGMPTrajectory(gmp, tau, y0, yg, pos_lim, vel_lim, accel_lim, true, false, Time, P_data, dP_data, ddP_data);
  // Data(Time, P_data, dP_data, ddP_data, "-", {0.85, 0.33, 0.1}, "opt-pos", true, true).write(fid, "5_");


  fid.write("pos_lim", pos_lim);
  fid.write("vel_lim", vel_lim);
  fid.write("accel_lim", accel_lim);
  fid.write("y0", y0);
  fid.write("yg", yg);
  fid.write("ygd", ygd);
  fid.write("tau", tau);
  fid.write("view", view_);

  fid.close();

  return 0;
}

// =================================================================
// =================================================================

void getGMPTrajectory(gmp_::GMP::Ptr gmp, double tau, const arma::vec &y0, const arma::vec &yg,
                arma::rowvec &Time, arma::mat &P_data, arma::mat &dP_data, arma::mat &ddP_data)
{
  Time.clear();
  P_data.clear();
  dP_data.clear();
  ddP_data.clear();

  arma::vec p = y0;
  arma::vec p_dot = arma::vec().zeros(p.size());
  arma::vec p_ddot = arma::vec().zeros(p.size());

  gmp->setY0(y0);
  gmp->setGoal(yg);

  double t = 0;
  double dt = 0.002;

  while (true)
  {
    double x = t/tau;
    double x_dot = 1/tau;
    double x_ddot = 0;

    if (x >= 1) x_dot = 0;

    arma::vec p_ref = gmp->getYd(x);
    arma::vec p_ref_dot = gmp->getYdDot(x, x_dot);
    arma::vec p_ref_ddot = gmp->getYdDDot(x, x_dot, x_ddot);

    arma::vec P = p_ref;
    arma::vec p_dot = p_ref_dot;
    arma::vec p_ddot = p_ref_ddot;

    Time = arma::join_horiz(Time, arma::vec({t}) );
    P_data = arma::join_horiz(P_data, p );
    dP_data = arma::join_horiz(dP_data, p_dot );
    ddP_data = arma::join_horiz(ddP_data, p_ddot );

    // p_ddot = p_ref_ddot + 30*(p_ref_dot - p_dot) + 100*(p_ref - p);

    t = t + dt;
    p = p + p_dot*dt;
    p_dot = p_dot + p_ddot*dt;

    if (x >= 1.0) break;
  }
}

void getOptGMPTrajectory(gmp_::GMP::Ptr gmp, double tau, const arma::vec &y0, const arma::vec &yg,
        const arma::mat &pos_lim, const arma::vec &vel_lim, const arma::vec &accel_lim, bool opt_pos, bool opt_vel,
        arma::rowvec &Time, arma::mat &P_data, arma::mat &dP_data, arma::mat &ddP_data)
{
  gmp_::GMP::Ptr gmp2( new gmp_::GMP() );
  gmp->deepCopy( gmp2.get() );


  gmp2->setScaleMethod( gmp_::TrajScale::Ptr( new gmp_::TrajScale_Prop(3) ) );
  gmp2->setY0(y0);
  gmp2->setGoal(yg);

  //Timer::tic();
 
  gmp_::GMP_Opt gmp_opt(gmp2.get());
  gmp_opt.setProblemOptions(opt_pos, opt_vel, false, 1, 1, 0.1);
  gmp_opt.setOptimizationOptions(2000, 0.035, true, true, false);
  gmp_opt.setMotionDuration(tau);

  int n_points = 100;
  arma::rowvec x_ineq = arma::linspace<arma::rowvec>(0,1, n_points);

  arma::mat pos_lb = arma::repmat(pos_lim.col(0), 1, n_points);
  arma::mat pos_ub = arma::repmat(pos_lim.col(1), 1, n_points);

  arma::mat vel_lb = arma::repmat(vel_lim(0)*arma::vec().ones(3), 1, n_points);
  arma::mat vel_ub = arma::repmat(vel_lim(1)*arma::vec().ones(3), 1, n_points);

  arma::mat accel_lb = arma::repmat(accel_lim(0)*arma::vec().ones(3), 1, n_points);
  arma::mat accel_ub = arma::repmat(accel_lim(1)*arma::vec().ones(3), 1, n_points);

  arma::rowvec x_eq = {0, 1};
  arma::mat pos_eq = arma::join_horiz(y0, yg);
  arma::mat vel_eq = arma::mat().zeros(3,2);
  arma::mat accel_eq = arma::mat().zeros(3,2);
/*

  // position constr
  gmp_opt.setPosBounds(pos_lim(:,1), pos_lim(:,2), n_points);
  gmp_opt.setPosConstr([],[],[], [0 1], [y0 yg]);
  // velocity constr
  gmp_opt.setVelBounds(vel_lim(1), vel_lim(2), n_points);
  gmp_opt.setVelConstr([], [], [], [0 1], zeros(3,2));
  // accel constr
  gmp_opt.setAccelBounds(accel_lim(1), accel_lim(2), n_points);
  gmp_opt.setAccelConstr([], [], [], [0 1], zeros(3,2));
*/
  gmp_opt.setPosConstr(x_ineq, pos_lb, pos_ub, x_eq, pos_eq);
  gmp_opt.setVelConstr(x_ineq, vel_lb, vel_ub, x_eq, vel_eq);
  gmp_opt.setAccelConstr(x_ineq, accel_lb, accel_ub, x_eq, accel_eq);

  // gmp_opt.optimize(100);
  Timer::tic();
  gmp_::GMP_Opt::SolutionStatus status = gmp_opt.optimize( arma::linspace<arma::rowvec>(0,1, 100) );
  std::cout << "Optimization finished: "; Timer::toc();
  if (status == gmp_::GMP_Opt::OPTIMAL) PRINT_INFO_MSG( gmp_opt.getExitMsg() + "\n" );
  else if (status == gmp_::GMP_Opt::SUBOPTIMAL) PRINT_WARNING_MSG( gmp_opt.getExitMsg() + "\n" );
  else /*if (status == gmp_::GMP_Opt::FAILED)*/ PRINT_ERROR_MSG( gmp_opt.getExitMsg() + "\n" );

  getGMPTrajectory(gmp2, tau, y0, yg, Time, P_data, dP_data, ddP_data);

}
