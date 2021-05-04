#include <cstdlib>
#include <memory>

#include <ros/ros.h>
#include <ros/package.h>

#include <armadillo>

#include <gmp_lib/gmp_lib.h>
#include <gmp_lib/io/file_io.h>

#include <project_name_/utils/utils.h>

using namespace as64_;

// ==================================
// ------------   MAIN  -------------
// ==================================

int main(int argc, char **argv)
{
  ros::init(argc, argv, "update_gmp_test_node");

  std::string gmp_filename = ros::package::getPath(PROJECT_NAME_) + "/src/optimization/data/gmp_pos.bin";

  // =============  Load GMP  =============
  gmp_::GMP_nDoF::Ptr gmp(new gmp_::GMP_nDoF(1, 2) );
  gmp_::GMP_nDoF_IO::read(gmp.get(), gmp_filename, "opt_");

  int n_dof = gmp->numOfDoFs();

  if (n_dof != 3) throw std::runtime_error("The number of DoFs must be equal to 3!\n");

  // =============  GMP update  =============
  double kt = 0.5;
  double ks = 2;
  arma::vec p0d = gmp->getYd(0); //Pd_data(:,1);
  arma::vec pgd = gmp->getYd(1); //Pd_data(:,end);
  arma::vec p0 = p0d;
  arma::vec pg = ks*(pgd-p0d) + p0;
  double T = 20/kt;
  Time = Timed/kt;
  x = Time / T;
  x_dot = 1/T;
  x_ddot = 0;

  gmp->setY0(p0);
  gmp->setGoal(pg);

  N = length(x);
  P_data = zeros(n_dof,N);
  dP_data = zeros(n_dof,N);
  ddP_data = zeros(n_dof,N);

  points = getPoint([], [], [], [], [], [], []);

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

  double t5 = 7 / kt;
  double x5 = t5/T;
  arma::vec p5_dot = gmp->getYdDot(x5,x_dot) + arma::vec({0.1, 0.05, -0.1});
  arma::vec p5_ddot = gmp->getYdDDot(x5,x_dot,x_ddot) + arma::vec({-0.03, 0.02, 0.025});

  gmp_::GMP_nDoF_Update::Ptr gmp_up(gmp);
  
  gmp_up.enableSigmawUpdate(true);
  gmp_up.setMsrNoiseVar(1e-4);
  for i=1:length(points), updateGMP(gmp_up, points(i)); end
  // updateGMP_all(gmp_up, points);

  // =============  Write results  =============
  // {
  //   gmp_::FileIO fid(results_filename, gmp_::FileIO::out | gmp_::FileIO::trunc);
  //   fid.write("Timed",Timed);
  //   fid.write("Pd_data",Pd_data);
  //   fid.write("dPd_data",dPd_data);
  //   fid.write("ddPd_data",ddPd_data);
  //   fid.write("Time",Time);
  //   fid.write("P_data",P_data);
  //   fid.write("dP_data",dP_data);
  //   fid.write("ddP_data",ddP_data);
  //   fid.write("spat_s",spat_s);
  //   fid.write("temp_s",temp_s);
  // }

  PRINT_CONFIRM_MSG("<====  Finished! ====>\n");

  return 0;
}
