#include <cstdlib>
#include <memory>

#include <ros/ros.h>
#include <ros/package.h>

#include <armadillo>

#include <gmp_lib/gmp_lib.h>
#include <gmp_lib/io/file_io.h>

#include <project_name_/utils/utils.h>

using namespace as64_;

std::string train_filename;
std::string results_filename;
std::string train_method = "LS";
int N_kernels = 25;
double kernels_std_scaling = 1.5;
std::string scale_type;
arma::vec wb_normal;
arma::vec spat_s;
double temp_s;

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
  gmp_::GMP_nDoF::Ptr gmp(new gmp_::GMP_nDoF(1, 2) );

  if (read_gmp_from_file) gmp_::GMP_nDoF_IO::read(gmp.get(), gmp_filename, "");
  else
  {
    // initialize and train GMP
    unsigned n_dof = Pd_data.n_rows;
    gmp.reset( new gmp_::GMP_nDoF(n_dof, N_kernels, kernels_std_scaling) );
    Timer::tic();
    arma::vec offline_train_mse;
    PRINT_INFO_MSG("GMP_nDoF training...\n");
    gmp->train(train_method, Timed/Timed.back(), Pd_data, &offline_train_mse);
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

  if (write_gmp_to_file) gmp_::GMP_nDoF_IO::write(gmp.get(), gmp_filename, "");

  // =============  GMP simulation  =============
  PRINT_INFO_MSG("GMP_nDoF simulation...\n");

  Timer::tic();

  int n_data = Timed.size();
  int i_end = n_data - 1;
  arma::vec P0d = Pd_data.col(0);
  arma::vec Pgd = Pd_data.col(i_end);
  arma::vec P0 = P0d;
  arma::vec Pg = spat_s%(Pgd - P0) + P0;
  double T = Timed(i_end) / temp_s;
  double dt = Ts;

  arma::rowvec Time;
  arma::mat P_data, dP_data, ddP_data;
  simulateGMP_nDoF(gmp, P0, Pg, T, dt, Time, P_data, dP_data, ddP_data);

  Timer::toc();

  // =============  Write results  =============
  {
    gmp_::FileIO fid(results_filename, gmp_::FileIO::out | gmp_::FileIO::trunc);
    fid.write("Timed",Timed);
    fid.write("Pd_data",Pd_data);
    fid.write("dPd_data",dPd_data);
    fid.write("ddPd_data",ddPd_data);
    fid.write("Time",Time);
    fid.write("P_data",P_data);
    fid.write("dP_data",dP_data);
    fid.write("ddP_data",ddP_data);
    fid.write("spat_s",spat_s);
    fid.write("temp_s",temp_s);
  }

  PRINT_CONFIRM_MSG("<====  Finished! ====>\n");

  return 0;
}

// ===============================================

void loadParams()
{
  std::string package_path = ros::package::getPath("project_name_") + "/";
  ros::NodeHandle nh("~");
  if (!nh.getParam("train_filename", train_filename)) throw std::runtime_error("Failed to load param \"train_filename\"...\n");
  train_filename = package_path + train_filename;
  if (!nh.getParam("results_filename", results_filename)) throw std::runtime_error("Failed to load param \"results_filename\"...\n");
  results_filename = package_path + results_filename;
  if (!nh.getParam("train_method", train_method)) throw std::runtime_error("Failed to load param \"train_method\"...\n");
  if (!nh.getParam("N_kernels", N_kernels)) throw std::runtime_error("Failed to load param \"N_kernels\"...\n");
  if (!nh.getParam("kernels_std_scaling", kernels_std_scaling)) throw std::runtime_error("Failed to load param \"kernels_std_scaling\"...\n");
  if (!nh.getParam("scale_type", scale_type)) throw std::runtime_error("Failed to load param \"scale_type\"...\n");
  std::vector<double> wb_normal_temp;
  if (!nh.getParam("wb_normal", wb_normal_temp)) throw std::runtime_error("Failed to load param \"wb_normal\"...\n");
  wb_normal = wb_normal_temp;

  std::vector<double> spat_s_;
  if (!nh.getParam("spat_s", spat_s_)) throw std::runtime_error("Failed to load param \"spat_s\"...\n");
  spat_s = spat_s_;
  if (!nh.getParam("temp_s", temp_s)) throw std::runtime_error("Failed to load param \"temp_s\"...\n");

  if (!nh.getParam("read_gmp_from_file", read_gmp_from_file)) throw std::runtime_error("Failed to load param \"read_gmp_from_file\"...\n");
  if (!nh.getParam("write_gmp_to_file", write_gmp_to_file)) throw std::runtime_error("Failed to load param \"write_gmp_to_file\"...\n");
  if ((read_gmp_from_file || write_gmp_to_file) && !nh.getParam("gmp_filename", gmp_filename))
    throw std::runtime_error("Failed to load param \"gmp_filename\"...\n");
  gmp_filename = package_path + gmp_filename;
}
