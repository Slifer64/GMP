#include <cstdlib>
#include <memory>

#include <ros/ros.h>
#include <ros/package.h>

#include <armadillo>

#include <gmp_lib/gmp_lib.h>
#include <gmp_lib/io/file_io.h>

#include <project_name_/utils/utils.h>

using namespace as64_;

typedef void (*sim_fun_ptr)(std::shared_ptr<gmp_::GMPo> &, const arma::vec &, const arma::vec &,
                            double , double , arma::rowvec &, arma::mat &, arma::mat &, arma::mat &);

void loadParams();

std::string path;

std::string train_data_file;
std::string sim_data_file;
std::string train_method;
unsigned N_kernels;
double D;
double K;
double kernels_std_scaling;
double ks;
double kt;
sim_fun_ptr simulateGMPo;

int main(int argc, char** argv)
{
  // ===========  Initialize the ROS node  ===============
  ros::init(argc, argv, "Test_orient_GMP_node");

  loadParams();

  // =============  Load train data  =============
  arma::rowvec Timed;
  arma::mat Qd_data, vRotd_data, dvRotd_data;
  gmp_::FileIO fid(train_filename, gmp_::FileIO::in);
  fid.read("Timed",Timed);
  fid.read("Qd_data",Qd_data);
  fid.read("vRotd_data",vRotd_data);
  fid.read("dvRotd_data",dvRotd_data);
  fid.close();

  double Ts = Timed(1) - Timed(0);

  // =============  Create/Train GMP  =============
  gmp_::GMPo::Ptr gmp(new gmp_::GMPo(2) );

  if (read_gmp_from_file) gmp_::GMPo_IO::read(gmp, gmp_filename, "");
  else
  {
    // initialize and train GMP
    unsigned n_dof = 3;
    gmp.reset( new gmp_::GMPo(N_kernels, kernels_std_scaling) );
    Timer::tic();
    arma::vec offline_train_mse;
    PRINT_INFO_MSG("GMPo training...\n");
    gmp->train(train_method, Timed/Timed.back(), Qd_data, &offline_train_mse);
    std::cerr << "offline_train_mse = \n" << offline_train_mse << "\n";
    Timer::toc();

    // set scaling type
    gmp_::TrajScale::Ptr traj_sc;
    if (scale_type.compare("prop") == 0) traj_sc.reset( new gmp_::TrajScale_Prop(n_dof) );
    else if (scale_type.compare("rot_min") == 0) traj_sc.reset( new gmp_::TrajScale_Rot_min() );
    else if (scale_type.compare("rot_wb") == 0)
    {
      traj_sc.reset( new gmp_::TrajScale_Rot_wb() );
      dynamic_cast<gmp_::TrajScale_Rot_wb *>(traj_sc.get())->setWorkBenchNormal( {0, 0, 1} );
    }
    else throw std::runtime_error("Unsupported scale type \"" + scale_type + "\"...\n");

    gmp->setScaleMethod(traj_sc);
  }

  if (write_gmp_to_file) gmp_::GMPo_IO::write(gmp, gmp_filename, "");

  /*
  // ===========  gmp update and simulation  ===============
  arma::vec Q0d = Qd_data.col(0);
  arma::vec Qgd = Qd_data.col(i_end);
  arma::vec Q0 = Q0d;
  arma::vec e0 = ks*gmp_::quatLog( gmp_::quatProd( Qgd, gmp_::quatInv(Q0d) ) );
  arma::vec Qg = gmp_::quatProd(gmp_::quatExp(e0), Q0);
  double T = kt*Timed(i_end);
  double dt = Timed(1) - Timed(0);

  arma::rowvec Time;
  arma::mat Q_data;
  arma::mat vRot_data;
  arma::mat dvRot_data;
  simulateGMPo(gmp, Q0, Qg, T, dt, Time, Q_data, vRot_data, dvRot_data);

  // ===========  write results  ===============
  io_::FileIO out(sim_data_file, io_::FileIO::out | io_::FileIO::trunc);
  out.write("Timed", Timed);
  out.write("Qd_data", Qd_data);
  out.write("vRotd_data", vRotd_data);
  out.write("dvRotd_data", dvRotd_data);
  out.write("Time", Time);
  out.write("Q_data", Q_data);
  out.write("vRot_data", vRot_data);
  out.write("dvRot_data", dvRot_data);
  out.write("ks", ks);
  out.write("kt", kt);
  */

  // ===========  Shutdown ROS node  ==================
  PRINT_CONFIRM_MSG("<====  Finished! ====>\n");

  return 0;
}

void loadParams()
{
  ros::NodeHandle nh_("~");

  // ===========  Read params  ===============
  path = ros::package::getPath("gmp_test") + "/matlab/data/";

  if (!nh_.getParam("train_data_file", train_data_file)) train_data_file = "gmp_train_data.bin";
  if (!nh_.getParam("sim_data_file", sim_data_file)) sim_data_file = "gmp_update_sim_data.bin";
  if (!nh_.getParam("train_method", train_method)) train_method = "LWR";
  int n_ker;
  if (!nh_.getParam("N_kernels", n_ker)) n_ker = 30;
  N_kernels = n_ker;
  if (!nh_.getParam("D", D)) D = 50;
  if (!nh_.getParam("K", K)) K = 250;
  if (!nh_.getParam("kernels_std_scaling", kernels_std_scaling)) kernels_std_scaling = 2;
  if (!nh_.getParam("ks", ks)) ks = 1.0;
  if (!nh_.getParam("kt", kt)) kt = 1.0;
  std::string sim_fun;
  if (!nh_.getParam("sim_fun", sim_fun)) sim_fun = "log";
  if (sim_fun.compare("log")==0) simulateGMPo = &simulateGMPo_in_log_space;
  else if (sim_fun.compare("quat")==0) simulateGMPo = &simulateGMPo_in_quat_space;
  else std::runtime_error("Unsupported simulation function: \"" + sim_fun + "\"...");

  train_data_file = path + train_data_file;
  sim_data_file = path + sim_data_file;
}
