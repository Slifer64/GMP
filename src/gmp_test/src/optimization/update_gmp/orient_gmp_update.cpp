#include <ros/ros.h>
#include <ros/package.h>

#include <armadillo>

#include <gmp_lib/gmp_lib.h>
#include <gmp_lib/math/quaternions.h>

#include <gmp_test/utils/utils.h>

using namespace as64_;

// ==================================
// ------------   MAIN  -------------
// ==================================

int main(int argc, char **argv)
{
  ros::init(argc, argv, "update_gmp_test_node");

  std::string in_filename = ros::package::getPath(PROJECT_NAME_) + "/src/optimization/update_gmp/data/gmp_orient.bin";
  std::string out_filename = ros::package::getPath(PROJECT_NAME_) + "/src/optimization/update_gmp/data/orient_results.bin";

  // =============  Load GMP  =============
  gmp_::GMPo::Ptr gmp_o( new gmp_::GMPo() );
  gmp_::read(gmp_o.get(), in_filename, "");
  
  double Ts = 0.01;
  arma::rowvec Timed = arma::linspace<arma::rowvec>(0,10, 10/Ts);
  arma::rowvec x = Timed/Timed.back();
  int n_data = x.size();

  // ==========    Construct new orientation trajectory

  // set demo and new init/target orientation
  arma::vec Q0d = gmp_o->getQd(0);
  arma::vec Qgd = gmp_o->getQd(1);

  arma::vec Q0 = Q0d;
  arma::vec Qg = gmp_::quatProd( gmp_::quatExp( arma::vec({0.5, 0.2, -0.4}) ), Qgd);

  // generate nominal trajectory for the new init/target orientation
  gmp_o->setQ0(Q0);
  gmp_o->setQg(Qg);
  arma::mat Qd_data(4,n_data);
  for (int j=0; j<n_data; j++) Qd_data.col(j) = gmp_o->getQd(x(j));

  arma::mat qd_data(3, n_data);

  arma::mat Qd_new(4, n_data);
  arma::mat qd_new_data(3, n_data);
  arma::mat vRotd_new_data(3, n_data);
  arma::mat dvRotd_new_data(3, n_data);

  for (int j=0; j<n_data; j++)
  {
    qd_data.col(j) = gmp_::GMPo::quat2q(Qd_data.col(j), Q0);
    qd_new_data.col(j) = qd_data.col(j) + 0.3*exp( -std::pow( (0.4-x(j))/0.1, 2) ) - 0.2*exp( -std::pow( (0.7-x(j))/0.06, 2) );
    Qd_new.col(j) = gmp_::GMPo::q2quat(qd_new_data.col(j), Q0);
  }

  for (int j=0; j<n_data-1; j++)
   vRotd_new_data.col(j) = gmp_::quatLog( gmp_::quatDiff(Qd_new.col(j+1),Qd_new.col(j)) )/Ts;

  for (int i=0;i<3;i++) dvRotd_new_data.row(i)= arma::join_horiz( arma::diff(vRotd_new_data.row(i)), arma::vec({0}) ) / Ts;

  // ==========    Update GMP

  double tau = Timed.back();
  double x_dot = 1/tau;
  double x_ddot = 0;

  gmp_::GMPo::Ptr gmp_o_new( new gmp_::GMPo() );
  gmp_o->deepCopy(gmp_o_new.get());
  gmp_::GMPo_Update gmp_up(gmp_o_new.get());
  gmp_up.enableSigmawUpdate(true);
  gmp_up.setMsrNoiseVar(1e-4);

  arma::uvec up_ind;
  double t_span = 0;
  for (int j=0; j<n_data; j++)
  {
    t_span += 1000*Ts;
    if ( arma::norm(qd_new_data.col(j)-qd_data.col(j) ) > 0.05  && t_span>300) // update every 300 ms
    {
      up_ind = arma::join_vert(up_ind, arma::uvec({j}));
      t_span = 0;
    }
  }
  x.elem(up_ind);
  arma::rowvec x_up = x.elem(up_ind).t();

  for (int j=0; j<up_ind.size(); j++)
  {
    int k = up_ind(j);
    gmp_up.updatePos(x_up(j), qd_new_data.col(k));
//    gmp_up.updateQuat(x_up(j), Qd_new.col(k));
//    gmp_up.updateRotVel(x_up(j), x_dot, vRotd_new_data.col(k), Qd_new.col(k));
//    gmp_up.updateRotAccel(x_up(j), x_dot, x_ddot, dvRotd_new_data.col(k), vRotd_new_data.col(k), Qd_new.col(k));
  }

  arma::mat q_data(3, n_data);
  arma::mat Q_data(4, n_data);
  arma::mat vRot_data(3, n_data);
  arma::mat dvRot_data(3, n_data);

  arma::mat q_new_data(3, n_data);
  arma::mat Q_new_data(4, n_data);
  arma::mat vRot_new_data(3, n_data);
  arma::mat dvRot_new_data(3, n_data);

  for (int j=0; j<n_data; j++)
  {
    q_data.col(j) = gmp_o->getYd(x(j));
    Q_data.col(j) = gmp_o->getQd(x(j));
    vRot_data.col(j) = gmp_o->getVd(x(j), x_dot);
    dvRot_data.col(j) = gmp_o->getVdDot(x(j), x_dot, x_ddot);

    q_new_data.col(j) = gmp_o_new->getYd(x(j));
    Q_new_data.col(j) = gmp_o_new->getQd(x(j));
    vRot_new_data.col(j) = gmp_o_new->getVd(x(j), x_dot);
    dvRot_new_data.col(j) = gmp_o_new->getVdDot(x(j), x_dot, x_ddot);
  }

  // =============  Write results  =============
  gmp_::FileIO fid(out_filename, gmp_::FileIO::out | gmp_::FileIO::trunc);
  fid.write("Timed", Timed);
  fid.write("q_data", q_data);
  fid.write("vRot_data", vRot_data);
  fid.write("dvRot_data", dvRot_data);
  fid.write("q_new_data", q_new_data);
  fid.write("vRot_new_data", vRot_new_data);
  fid.write("dvRot_new_data", dvRot_new_data);
  fid.write("qd_new_data", qd_new_data);
  fid.write("vRotd_new_data", vRotd_new_data);
  fid.write("dvRotd_new_data", dvRotd_new_data);
  fid.write("up_ind", up_ind);
  fid.close();

  PRINT_CONFIRM_MSG("<====  Finished! ====>\n");

  return 0;
}