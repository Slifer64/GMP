#include <ros/ros.h>
#include <ros/package.h>

#include <armadillo>

#include <gmp_lib/gmp_lib.h>

#include <gmp_test/utils/utils.h>

using namespace as64_;

// ==================================
// ------------   MAIN  -------------
// ==================================

int main(int argc, char **argv)
{
  ros::init(argc, argv, "update_gmp_test_node");

  std::string in_filename = ros::package::getPath(PROJECT_NAME_) + "/src/optimization/update_gmp/gmp_pos.bin";
  std::string out_filename = ros::package::getPath(PROJECT_NAME_) + "/src/optimization/update_gmp/results.bin";

  // =============  Load GMP  =============
  gmp_::GMP::Ptr gmp( new gmp_::GMP() );
  gmp_::read(gmp.get(), in_filename, "up_");

  int n_dof = gmp->numOfDoFs();

  if (n_dof != 3) throw_error("The number of DoFs must be equal to 3!\n");

  // =============  GMP update  =============
  double kt = 0.5;
  double ks = 2;

  double T = 10/kt;
  double x_dot = 1/T;
  double x_ddot = 0;

  arma::vec p0d = gmp->getYd(0); //Pd_data(:,1);
  arma::vec pgd = gmp->getYd(1); //Pd_data(:,end);
  arma::vec p0 = p0d;
  arma::vec pg = ks*(pgd-p0d) + p0;
  gmp->setY0(p0);
  gmp->setGoal(pg);

  double t1 = 1.5 / kt;
  double x1 = t1/T;
  arma::vec p1 = gmp->getYd(x1) + arma::vec({-0.1, 0.08, -0.12});

  double t2 = 3 / kt;
  double x2 = t2/T;
  arma::vec p2_dot = gmp->getYdDot(x2,x_dot) + arma::vec({0.05, -0.08, 0.1});

  double t3 = 4.25 / kt;
  double x3 = t3/T;
  arma::vec p3_ddot = gmp->getYdDDot(x3,x_dot,x_ddot) + arma::vec({-0.02, 0.03, 0.05});

  double t4 = 5.5 / kt;
  double x4 = t4/T;
  arma::vec p4 = gmp->getYd(x4) + arma::vec({0.2, -0.15, -0.13});
  arma::vec p4_ddot = gmp->getYdDDot(x4,x_dot,x_ddot) + arma::vec({0.02, -0.04, -0.025});
  std::vector<gmp_::Phase> s4 = { gmp_::Phase(x4,x_dot,x_ddot), gmp_::Phase(x4,x_dot,x_ddot) };
  std::vector<gmp_::UPDATE_TYPE> type4 = { gmp_::UPDATE_TYPE::POS, gmp_::UPDATE_TYPE::ACCEL };
  arma::mat Z4 = arma::join_horiz(p4, p4_ddot);

  double t5 = 7 / kt;
  double x5 = t5/T;
  arma::vec p5_dot = gmp->getYdDot(x5,x_dot) + arma::vec({0.1, 0.05, -0.1});
  arma::vec p5_ddot = gmp->getYdDDot(x5,x_dot,x_ddot) + arma::vec({-0.03, 0.02, 0.025});
  std::vector<gmp_::Phase> s5 = { gmp_::Phase(x5,x_dot,x_ddot), gmp_::Phase(x5,x_dot,x_ddot) };
  std::vector<gmp_::UPDATE_TYPE> type5 = { gmp_::UPDATE_TYPE::VEL, gmp_::UPDATE_TYPE::ACCEL };
  arma::mat Z5 = arma::join_horiz(p5_dot, p5_ddot);

  // make a deep copy of the original GMP
  gmp_::GMP::Ptr gmp0( new gmp_::GMP() );
  gmp->deepCopy(gmp0.get());

  gmp_::GMP_Update gmp_up(gmp.get());
  //gmp_up.initSigmaWfromMsr(arma::linspace<arma::rowvec>(0,1,100));
  gmp_up.enableSigmawUpdate(true);
  gmp_up.setMsrNoiseVar(1e-4);

  gmp_up.updatePos(x1, p1);
  gmp_up.updateVel(x2, x_dot, p2_dot);
  gmp_up.updateAccel(x3, x_dot, x_ddot, p3_ddot);
  gmp_up.updateWeights(s4, Z4, type4);
  gmp_up.updateWeights(s5, Z5, type5);

  // =============  Generate original and new trajectories  =============
  double Ts = 0.01;
  arma::rowvec Time = arma::linspace<arma::rowvec>(0,T, std::round(T/Ts));
  arma::rowvec x = Time/Time.back();
  int n_data = x.size();

  arma::mat P_data(3, n_data);
  arma::mat dP_data(3, n_data);
  arma::mat ddP_data(3, n_data);

  arma::mat P_new_data(3, n_data);
  arma::mat dP_new_data(3, n_data);
  arma::mat ddP_new_data(3, n_data);
  for (int j=0; j<n_data; j++)
  {
    P_data.col(j) = gmp0->getYd(x(j));
    dP_data.col(j) = gmp0->getYdDot(x(j), x_dot);
    ddP_data.col(j) = gmp0->getYdDDot(x(j), x_dot, x_ddot);

    P_new_data.col(j) = gmp->getYd(x(j));
    dP_new_data.col(j) = gmp->getYdDot(x(j), x_dot);
    ddP_new_data.col(j) = gmp->getYdDDot(x(j), x_dot, x_ddot);
  }

  // =============  Write results  =============
  gmp_::FileIO fid(out_filename, gmp_::FileIO::out | gmp_::FileIO::trunc);
  fid.write("Time", Time);
  fid.write("P_data", P_data);
  fid.write("dP_data", dP_data);
  fid.write("ddP_data", ddP_data);
  fid.write("P_new_data", P_new_data);
  fid.write("dP_new_data", dP_new_data);
  fid.write("ddP_new_data", ddP_new_data);
  fid.write("p1", p1);
  fid.write("p2_dot", p2_dot);
  fid.write("p3_ddot", p3_ddot);
  fid.write("p4", p4);
  fid.write("p4_ddot", p4_ddot);
  fid.write("p5_dot", p5_dot);
  fid.write("p5_ddot", p5_ddot);
  fid.close();

  PRINT_CONFIRM_MSG("<====  Finished! ====>\n");

  return 0;
}
