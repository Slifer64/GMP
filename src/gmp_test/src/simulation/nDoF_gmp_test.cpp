#include <cstdlib>
#include <memory>

#include <ros/ros.h>
#include <ros/package.h>

#include <armadillo>

#include <gmp_lib/gmp_lib.h>

#include <gmp_test/utils/utils.h>

using namespace as64_;

std::string train_filename;
std::string results_filename;
std::string train_method;
int N_kernels;
double kernels_std_scaling;
std::string scale_type;
arma::vec wb_normal;
arma::vec Pg;
double T;

bool read_gmp_from_file;
bool write_gmp_to_file;
std::string gmp_filename;

void loadParams();

// ==================================
// ------------   MAIN  -------------
// ==================================

int main(int argc, char **argv)
{
  ros::init(argc, argv, "nDoF_gmp_test_node");

  // =============  Load params  =============
  loadParams();

  // =============  Load train data  =============
  arma::rowvec Timed;
  arma::mat Pd_data, dPd_data, ddPd_data;
  gmp_::FileIO fid(train_filename, gmp_::FileIO::in);
  fid.read("Timed",Timed);
  fid.read("Pd_data",Pd_data);
  fid.read("dPd_data",dPd_data);
  fid.read("ddPd_data",ddPd_data);
  fid.close();

  double Ts = Timed(1) - Timed(0);

  // =============  Create/Train GMP  =============
  gmp_::GMP::Ptr gmp( new gmp_::GMP() );

  unsigned n_dof = Pd_data.n_rows;
  if (read_gmp_from_file)
  {
    gmp_::read(gmp.get(), gmp_filename, "");
    PRINT_INFO_MSG("Loaded GMP from file!\n");
  }
  else
  {
    // initialize and train GMP
    gmp.reset( new gmp_::GMP(n_dof, N_kernels, kernels_std_scaling) );
    Timer::tic();
    arma::vec offline_train_mse;
    PRINT_INFO_MSG("GMP training...\n");
    gmp->train(train_method, Timed/Timed.back(), Pd_data, &offline_train_mse);
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

  gmp->setScaleMethod(traj_sc);

  if (write_gmp_to_file)
  {
    gmp_::write(gmp.get(), gmp_filename, "");
    PRINT_INFO_MSG("Wrote GMP to file!\n");
  }

  // gmp->autoRetrain(50, 2, 200, "LWR");

  // =============  GMP simulation  =============
  PRINT_INFO_MSG("GMP simulation...\n");

  Timer::tic();

  int n_data = Timed.size();
  int i_end = n_data - 1;
  arma::vec P0d = Pd_data.col(0);
  arma::vec Pgd = Pd_data.col(i_end);
  arma::vec P0 = P0d;
  // arma::vec Pg = spat_s%(Pgd - P0) + P0;
  // double T = Timed(i_end) / temp_s;
  double dt = Ts;

  arma::rowvec Time;
  arma::mat P_data, dP_data, ddP_data;
  simulateGMP(gmp, P0, Pg, T, dt, Time, P_data, dP_data, ddP_data);

  Timer::toc();

  // This is the groundtruth trajectory that should be produced
  arma::mat Ks = gmp->getScaling();
  double temp_s = Timed.back() / T; // temporal scaling 
  // arma::rowvec Timed2 = Timed / temp_s;
  // arma::mat Pd2_data = Ks*( Pd_data-P0d ) + P0;
  // arma::mat dPd2_data = Ks*dPd_data*temp_s;
  // arma::mat ddPd2_data = Ks*ddPd_data*std::pow(temp_s,2);

  // =============  Write results  =============
  {
    gmp_::FileIO fid(results_filename, gmp_::FileIO::out | gmp_::FileIO::trunc);
    // write the demo data
    fid.write("Timed",Timed);
    fid.write("Pd_data",Pd_data);
    fid.write("dPd_data",dPd_data);
    fid.write("ddPd_data",ddPd_data);
    // write the DMP data
    fid.write("Time",Time);
    fid.write("P_data",P_data);
    fid.write("dP_data",dP_data);
    fid.write("ddP_data",ddP_data);
    // write the spatial scaling matrix and the temporal scaling
    fid.write("Ks",Ks);
    fid.write("temp_s",temp_s);
  }

  PRINT_CONFIRM_MSG("<====  Finished! ====>\n");

  return 0;
}

// ===============================================

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

  std::vector<double> Pg_;
  if (!nh.getParam("Pg", Pg_)) throw_error("Failed to load param \"Pg\"...\n");
  Pg = Pg_;
  if (!nh.getParam("T", T)) throw_error("Failed to load param \"T\"...\n");

  if (!nh.getParam("read_gmp_from_file", read_gmp_from_file)) throw_error("Failed to load param \"read_gmp_from_file\"...\n");
  if (!nh.getParam("write_gmp_to_file", write_gmp_to_file)) throw_error("Failed to load param \"write_gmp_to_file\"...\n");
  if ((read_gmp_from_file || write_gmp_to_file) && !nh.getParam("gmp_filename", gmp_filename))
    throw_error("Failed to load param \"gmp_filename\"...\n");
  gmp_filename = package_path + gmp_filename;
}
