#include <cstdlib>
#include <memory>

#include <ros/ros.h>
#include <ros/package.h>

#include <armadillo>

#include <gmp_lib/gmp_lib.h>

#include <gmp_test/utils/utils.h>

using namespace as64_;

// ##############################################################

typedef void (*sim_fun_ptr)(std::shared_ptr<gmp_::GMPo> &, const arma::vec &, const arma::vec &,
                            double , double , arma::rowvec &, arma::mat &, arma::mat &, arma::mat &);

void loadParams();

// ##############################################################

std::string train_filename;
std::string results_filename;
std::string train_method = "LS";
int N_kernels = 25;
double kernels_std_scaling = 1.5;
std::string scale_type;
arma::vec wb_normal;
arma::vec spat_s;
double temp_s;
sim_fun_ptr simulateGMPo;

bool read_gmp_from_file;
bool write_gmp_to_file;
std::string gmp_filename;

// ##############################################################

int main(int argc, char** argv)
{
  // ===========  Initialize the ROS node  ===============
  ros::init(argc, argv, "orient_gmp_test_node");

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
  unsigned n_dof = 3;
  gmp_::GMPo::Ptr gmp_o(new gmp_::GMPo(2) );

  if (read_gmp_from_file) gmp_::read(gmp_o.get(), gmp_filename, "");
  else
  {
    // initialize and train GMP
    gmp_o.reset( new gmp_::GMPo(N_kernels, kernels_std_scaling) );
    Timer::tic();
    arma::vec offline_train_mse;
    PRINT_INFO_MSG("GMPo training...\n");
    gmp_o->train(train_method, Timed/Timed.back(), Qd_data, &offline_train_mse);
    std::cerr << "offline_train_mse = \n" << offline_train_mse << "\n";
    Timer::toc();
  }

  // set scaling type
  gmp_::TrajScale::Ptr traj_sc;
  if (scale_type.compare("prop") == 0) traj_sc.reset( new gmp_::TrajScale_Prop(n_dof) );
  else if (scale_type.compare("rot_min") == 0) traj_sc.reset( new gmp_::TrajScale_Rot_min() );
  else if (scale_type.compare("rot_wb") == 0)
  {
    traj_sc.reset( new gmp_::TrajScale_Rot_wb() );
    dynamic_cast<gmp_::TrajScale_Rot_wb *>(traj_sc.get())->setWorkBenchNormal( {0, 0, 1} );
  }
  else throw_error("Unsupported scale type \"" + scale_type + "\"...\n");

  gmp_o->setScaleMethod(traj_sc);

  if (write_gmp_to_file) gmp_::write(gmp_o.get(), gmp_filename, "");

  // =============  GMP simulation  =============
  PRINT_INFO_MSG("GMP_orient simulation...\n");

  Timer::tic();

  int i_end = Timed.size() - 1;
  arma::vec Qd0 = Qd_data.col(0);
  arma::vec Q0 = Qd0; // math_.quatProd(rotm2quat(rotz(57))',Qd0);
  arma::vec Qgd = Qd_data.col(i_end);

  arma::mat ks = arma::diagmat(spat_s);
  double kt = temp_s;
  arma::vec e0 = ks*gmp_::quatLog(gmp_::quatDiff(Qgd,Qd0));
  arma::vec Qg = gmp_::quatProd(gmp_::quatExp(e0), Q0);
  double T = Timed(i_end) / kt;
  double dt = Ts;

  arma::rowvec Time;
  arma::mat Q_data, vRot_data, dvRot_data;
  simulateGMPo(gmp_o, Q0, Qg, T, dt, Time, Q_data, vRot_data, dvRot_data);

  Timer::toc();

  // =============  Write results  =============
  {
    gmp_::FileIO fid(results_filename, gmp_::FileIO::out | gmp_::FileIO::trunc);
    fid.write("Timed", Timed);
    fid.write("Qd_data", Qd_data);
    fid.write("vRotd_data", vRotd_data);
    fid.write("dvRotd_data", dvRotd_data);
    fid.write("Time", Time);
    fid.write("Q_data", Q_data);
    fid.write("vRot_data", vRot_data);
    fid.write("dvRot_data", dvRot_data);
    arma::mat T_sc = gmp_o->getScaling();
    fid.write("T_sc", T_sc);
    fid.write("temp_s", kt);
  }

  // ===========  Shutdown ROS node  ==================
  PRINT_CONFIRM_MSG("<====  Finished! ====>\n");

  return 0;
}

// ##############################################################

void loadParams()
{
  std::string package_path = ros::package::getPath(PROJECT_NAME_) + "/";
  ros::NodeHandle nh("~");
  if (!nh.getParam("train_filename", train_filename)) throw_error("Failed to load param \"train_filename\"...\n");
  train_filename = package_path + train_filename;
  if (!nh.getParam("results_filename", results_filename)) throw_error("Failed to load param \"results_filename\"...\n");
  results_filename = package_path + results_filename;
  if (!nh.getParam("train_method", train_method)) throw_error("Failed to load param \"train_method\"...\n");
  if (!nh.getParam("N_kernels", N_kernels)) throw_error("Failed to load param \"N_kernels\"...\n");
  if (!nh.getParam("kernels_std_scaling", kernels_std_scaling)) throw_error("Failed to load param \"kernels_std_scaling\"...\n");
  if (!nh.getParam("scale_type", scale_type)) throw_error("Failed to load param \"scale_type\"...\n");
  std::vector<double> wb_normal_temp;
  if (!nh.getParam("wb_normal", wb_normal_temp)) throw_error("Failed to load param \"wb_normal\"...\n");
  wb_normal = wb_normal_temp;

  std::vector<double> spat_s_;
  if (!nh.getParam("spat_s", spat_s_)) throw_error("Failed to load param \"spat_s\"...\n");
  spat_s = spat_s_;
  if (!nh.getParam("temp_s", temp_s)) throw_error("Failed to load param \"temp_s\"...\n");

  std::string orient_sim_fun;
  if (!nh.getParam("orient_sim_fun", orient_sim_fun)) throw_error("Failed to load param \"orient_sim_fun\"...\n");
  if (orient_sim_fun.compare("log") == 0) simulateGMPo = &simulateGMPo_in_log_space;
  else if (orient_sim_fun.compare("quat") == 0) simulateGMPo = &simulateGMPo_in_quat_space;
  else if (orient_sim_fun.compare("Cart") == 0) simulateGMPo = &simulateGMPo_in_Cart_space;
  else throw_error("Unsupported orient-simulation function \"" + orient_sim_fun + "\"...\n");

  if (!nh.getParam("read_gmp_from_file", read_gmp_from_file)) throw_error("Failed to load param \"read_gmp_from_file\"...\n");
  if (!nh.getParam("write_gmp_to_file", write_gmp_to_file)) throw_error("Failed to load param \"write_gmp_to_file\"...\n");
  if ((read_gmp_from_file || write_gmp_to_file) && !nh.getParam("gmp_filename", gmp_filename))
    throw_error("Failed to load param \"gmp_filename\"...\n");
  gmp_filename = package_path + gmp_filename;

}
