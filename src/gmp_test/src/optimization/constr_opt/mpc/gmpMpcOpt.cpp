#include <iostream>

#include <gmp_lib/GMP/GMP_MPC.h>
#include <gmp_lib/io/file_io.h>

#include <ros/package.h>

using namespace as64_;

void PRINT_INFO_MSG(const std::string &msg, std::ostream &out=std::cout)
{
  out << "\033[1m" << "\033[34m" << "[INFO]: " << msg << "\033[0m";
}

void PRINT_CONFIRM_MSG(const std::string &msg, std::ostream &out=std::cout)
{
  out << "\033[1m" << "\033[32m" << "[INFO]: " << msg << "\033[0m";
}

void PRINT_WARNING_MSG(const std::string &msg, std::ostream &out=std::cerr)
{
  out << "\033[1m" << "\033[33m" << "[WARNING]: " << msg << "\033[0m";
}

void PRINT_ERROR_MSG(const std::string &msg, std::ostream &out=std::cerr)
{
  out << "\033[1m" << "\033[31m" << "[ERROR]: " << msg << "\033[0m";
}

int main(int argc, char **argv)
{
    gmp_::GMP::Ptr gmp(new gmp_::GMP());
    gmp_::GMP_MPC gmp_mpc( gmp.get(), 10, 0.02, 30, 1.5, {0,0,0}, {1e6,100,20});

    return 0;
}

void gmpMpcOpt(const gmp_::GMP *gmp0, double tau, const arma::vec &y0, const arma::vec &yg,
               const arma::mat &pos_lim, const arma::mat &vel_lim, const arma::mat &accel_lim)
{
  std::string filename = ros::package::getPath("gmp_test") + "/optimization/constr_opt/mpc/data/input.bin";
  gmp_::FileIO fid(filename, gmp_::FileIO::in);

  double dt;
  fid.read("dt",dt);

  double tau;
  fid.read("tau",tau);

  unsigned n_dof = y0.size();

  arma::rowvec Time;
  arma::mat P_data;
  arma::mat dP_data;
  arma::mat ddP_data;
  arma::mat slack_data;

  arma::vec O_ndof = arma::vec().zeros(n_dof);

  double t = 0;
  double s = t/tau;
  double s_dot = 1/tau;
  double s_ddot = 0;
  arma::vec y = y0;
  arma::vec y_dot = O_ndof;
  arma::vec y_ddot = O_ndof;

  bool pos_slack = 1;
  bool vel_slack = 1;
  bool accel_slack = 1;
  std::vector<bool> slack_flags = {pos_slack, vel_slack, accel_slack};
  std::vector<double> slack_gains = {1e6, 100, 20};

  // --------  GMP - MPC  --------
  gmp_::GMP_MPC gmp_mpc(gmp.get(), 10, 0.02, 30, 1.5, slack_flags, slack_gains);

  gmp_mpc.setPosLimits(pos_lim.col(0), pos_lim.col(1));
  gmp_mpc.setVelLimits(vel_lim.col(0), vel_lim.col(1));
  gmp_mpc.setAccelLimits(accel_lim.col(0), accel_lim.col(1));

  // gmp_mpc.setPosSlackLimit(1e-3);
  // gmp_mpc.setVelSlackLimit(0.05);
  // gmp_mpc.setAccelSlackLimit(0.1);

  gmp_mpc.setPosSlackLimit(0);
  gmp_mpc.setVelSlackLimit(0);
  gmp_mpc.setAccelSlackLimit(0);

  gmp_mpc.setInitialState(y, y_dot, y_ddot, s, s_dot, s_ddot);
  gmp_mpc.setFinalState(yg, O_ndof, O_ndof, 1, 0, 0);

  gmp->setScaleMethod( gmp_::TrajScale::Ptr( new gmp_::TrajScale_Prop(n_dof) ) );
  gmp->setY0(y0);
  gmp->setGoal(yg);

  arma::wall_clock timer;
  timer.tic();
  // --------  Simulation loop  --------
  while (true)
  {
    // --------  Stopping criteria  --------
    if (s > 1.0) break;

    //text_prog.update(100*t/tau);

    if (s >= 1)
    {
      s = 1;
      s_dot = 0;
      s_ddot = 0;
    }

    gmp_::GMP_MPC::Solution sol = gmp_mpc.solve(s, s_dot, s_ddot);
    arma::vec y = sol.y;
    arma::vec y_dot = sol.y_dot;
    arma::vec y_ddot = sol.y_ddot;
    arma::vec slack_var = sol.slack_var;

    // gmp_mpc.setInitialState(y, y_dot, y_ddot, s, s_dot, s_ddot);

    if (sol.exit_flag > 0) PRINT_WARNING_MSG(sol.exit_msg);
    else if (sol.exit_flag < 0)
    {
      PRINT_ERROR_MSG(sol.exit_msg);
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
}

