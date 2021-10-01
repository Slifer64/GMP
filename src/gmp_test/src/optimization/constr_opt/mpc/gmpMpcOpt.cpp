#include <iostream>
#include <algorithm>

#include <gmp_lib/GMP/GMP_MPC.h>
#include <gmp_lib/io/gmp_io.h>

#include <ros/package.h>

using namespace as64_;

void PRINT_INFO_MSG(const std::string &msg, std::ostream &out=std::cout);
void PRINT_CONFIRM_MSG(const std::string &msg, std::ostream &out=std::cout);
void PRINT_WARNING_MSG(const std::string &msg, std::ostream &out=std::cerr);
void PRINT_ERROR_MSG(const std::string &msg, std::ostream &out=std::cerr);

int main(int argc, char **argv)
{
  // data path
  std::string path = ros::package::getPath("gmp_test") + "/src/optimization/constr_opt/mpc/data/";

  // read input data
  arma::vec y0;
  arma::vec yg;
  double tau;
  arma::mat pos_lim;
  arma::mat vel_lim;
  arma::mat accel_lim;
  gmp_::GMP::Ptr gmp(new gmp_::GMP());
  {
    gmp_::FileIO fid(path + "gmp_mpc_opt_in.bin", gmp_::FileIO::in);
    fid.read("y0", y0);
    fid.read("yg", yg);
    fid.read("tau", tau);
    fid.read("pos_lim", pos_lim);
    fid.read("vel_lim", vel_lim);
    fid.read("accel_lim", accel_lim);
    gmp_::read(gmp.get(), fid, "gmp_");
  }

  unsigned n_dof = y0.size();
  arma::vec ones_ndof = arma::vec().ones(n_dof);
  arma::vec O_ndof = arma::vec().zeros(n_dof);

  arma::rowvec Time;
  arma::mat P_data;
  arma::mat dP_data;
  arma::mat ddP_data;
  arma::mat slack_data;

  double t = 0;
  double dt = 0.002;
  double s = t/tau;
  double s_dot = 1/tau;
  double s_ddot = 0;
  arma::vec y = y0;
  arma::vec y_dot = O_ndof;
  arma::vec y_ddot = O_ndof;

  // --------  GMP - MPC  --------
  int N_horizon = 12;
  double pred_time_step = 0.02;
  int N_kernels = 30;
  double kernels_std_scaling = 1.5;

  bool pos_slack = 1;
  bool vel_slack = 1;
  bool accel_slack = 1;
  std::vector<bool> slack_flags = {pos_slack, vel_slack, accel_slack};
  std::vector<double> slack_gains = {1e6, 100, 20};

  gmp_::GMP_MPC gmp_mpc(gmp.get(), N_horizon, pred_time_step, N_kernels, kernels_std_scaling, slack_flags, slack_gains);

  gmp_mpc.settings.time_limit = 1.5e-3;
//  gmp_mpc.settings.max_iter = 4000;
  gmp_mpc.settings.abs_tol = 1e-3;
  gmp_mpc.settings.rel_tol = 1e-4;

  gmp_mpc.setPosLimits(pos_lim.col(0), pos_lim.col(1));
  gmp_mpc.setVelLimits(vel_lim.col(0), vel_lim.col(1));
  gmp_mpc.setAccelLimits(accel_lim.col(0), accel_lim.col(1));

  gmp_mpc.setPosSlackLimit(1e-3);
  gmp_mpc.setVelSlackLimit(0.02);
  gmp_mpc.setAccelSlackLimit(0.2);

  gmp_mpc.setInitialState(y, y_dot, y_ddot, s, s_dot, s_ddot);
  gmp_mpc.setFinalState(yg, O_ndof, O_ndof, 1, 0, 0);

  gmp->setScaleMethod( gmp_::TrajScale::Ptr( new gmp_::TrajScale_Prop(n_dof) ) );
  gmp->setY0(y0);
  gmp->setGoal(yg);

  arma::rowvec elaps_t_data;

  arma::wall_clock timer, timer2;
  timer.tic();
  // --------  Simulation loop  --------
  while (true)
  {
    // --------  Stopping criteria  --------
    if (s > 1.0) break;

    //text_prog.update(100*t/tau);

    // std::cout << 100*t/tau << "\n";

    if (s >= 1)
    {
      s = 1;
      s_dot = 0;
      s_ddot = 0;
    }

    if (s > 0.9)
    {
      gmp_mpc.settings.abs_tol = 1e-4;
      gmp_mpc.settings.rel_tol = 1e-5;
    }


    timer2.tic();
    gmp_::GMP_MPC::Solution sol = gmp_mpc.solve(s, s_dot, s_ddot);
    arma::vec y = sol.y;
    arma::vec y_dot = sol.y_dot;
    arma::vec y_ddot = sol.y_ddot;
    arma::vec slack_var = sol.slack_var;
    elaps_t_data = arma::join_horiz( elaps_t_data, arma::vec({timer2.toc()*1000}) );

    // gmp_mpc.setInitialState(y, y_dot, y_ddot, s, s_dot, s_ddot);

    if (sol.exit_flag > 0) PRINT_WARNING_MSG(sol.exit_msg + "\n");
    else if (sol.exit_flag < 0)
    {
      PRINT_ERROR_MSG(sol.exit_msg + "\n");
      break;
    }

    // --------  Log data  --------
    Time = arma::join_horiz(Time, arma::mat({t}));
    P_data = arma::join_horiz(P_data, y);
    dP_data = arma::join_horiz(dP_data, y_dot);
    ddP_data = arma::join_horiz(ddP_data, y_ddot);
    slack_data = arma::join_horiz(slack_data, slack_var);

    // --------  Numerical integration  --------
    t = t + dt;
    s = s + s_dot*dt;
    s_dot = s_dot + s_ddot*dt;
  }

  //text_prog.update(100);

  double elaps_t_ms = timer.toc()*1000;
  PRINT_INFO_MSG("===> GMP-MPC optimization finished! Elaps time: " + std::to_string(elaps_t_ms) + " ms\n");

  double mean_elaps_t = arma::mean(elaps_t_data);
  double std_elaps_t = arma::stddev(elaps_t_data);
  double max_elaps_t = arma::max(elaps_t_data);
  double min_elaps_t = arma::min(elaps_t_data);

  std::cout << "======= Elapsed time (ms) ======\n";
  std::cout << "std_range: [" << std::max(0.,mean_elaps_t-std_elaps_t) << " -  " << mean_elaps_t + std_elaps_t <<"]\n";
  std::cout << "mean     : " << mean_elaps_t << " +/- " << std_elaps_t <<"\n";
  std::cout << "min      : " << min_elaps_t << "\n";
  std::cout << "max      : " << max_elaps_t << "\n";
  std::cout << "==============================\n";

  // ----------- write results ----------------
  {
      gmp_::FileIO fid(path + "gmp_mpc_opt_out.bin", gmp_::FileIO::out | gmp_::FileIO::trunc);
      fid.write("Time", Time);
      fid.write("P_data", P_data);
      fid.write("dP_data", dP_data);
      fid.write("ddP_data", ddP_data);
      fid.write("slack_data", slack_data);
      fid.write("pos_slack", pos_slack);
      fid.write("vel_slack", vel_slack);
      fid.write("accel_slack", accel_slack);
  }


  return 0;
}

void PRINT_INFO_MSG(const std::string &msg, std::ostream &out)
{
  out << "\033[1m" << "\033[34m" << "[INFO]: " << msg << "\033[0m";
}

void PRINT_CONFIRM_MSG(const std::string &msg, std::ostream &out)
{
  out << "\033[1m" << "\033[32m" << "[INFO]: " << msg << "\033[0m";
}

void PRINT_WARNING_MSG(const std::string &msg, std::ostream &out)
{
  out << "\033[1m" << "\033[33m" << "[WARNING]: " << msg << "\033[0m";
}

void PRINT_ERROR_MSG(const std::string &msg, std::ostream &out)
{
  out << "\033[1m" << "\033[31m" << "[ERROR]: " << msg << "\033[0m";
}

