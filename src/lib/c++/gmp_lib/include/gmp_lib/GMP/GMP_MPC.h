#ifndef GMP_LIB_GMP_MPC_H
#define GMP_LIB_GMP_MPC_H

// N-DoF GMP-MPC Optimization class
//

#include <gmp_lib/GMP/GMP.h>
#include <list>

#include <armadillo>

namespace as64_
{

namespace gmp_
{

class GMP_MPC
{

public:

  typedef std::shared_ptr<GMP_MPC> Ptr;

  struct
  {
    unsigned max_iter;
    double time_limit;
    double abs_tol;
    double rel_tol;
  } settings;

  struct Solution
  {
    arma::vec y;
    arma::vec y_dot;
    arma::vec y_ddot;

    arma::vec pos_slack;
    arma::vec vel_slack;
    arma::vec accel_slack;

    int exit_flag;
    std::string exit_msg;
  };
  
  GMP_MPC(const GMP *gmp, unsigned N_horizon, double pred_time_step, unsigned N_kernels, double kernel_std_scaling, const std::vector<double> &slack_gains);

  void setInitialState(const arma::vec &y0, const arma::vec &y0_dot, const arma::vec &y0_ddot, double s, double s_dot, double s_ddot);
  
  void setFinalState(const arma::vec &yf, const arma::vec &yf_dot, const arma::vec &yf_ddot, double s, double s_dot, double s_ddot);
  void setFinalState(const arma::vec &yf, const arma::vec &yf_dot, const arma::vec &yf_ddot, double s, double s_dot, double s_ddot, const arma::vec &err_tol);

  void setObjCostGains(double pos_gain, double vel_gain);

  void setPosLimits(const arma::vec &lb, const arma::vec &ub);
  
  void setVelLimits(const arma::vec &lb, const arma::vec &ub);
  
  void setAccelLimits(const arma::vec &lb, const arma::vec &ub);
  
  void addViaPoint(double s, const arma::vec &y, double err_tol = 1e-6);

  void setPosSlackLimit(double s_lim);
  
  void setVelSlackLimit(double s_lim);
  
  void setAccelSlackLimit(double s_lim);

  GMP_MPC::Solution solve(double s, double s_dot, double s_ddot);

  arma::vec getYd(double s) const;

  arma::vec getYdDot(double s, double s_dot) const;

  arma::vec getYdDDot(double s, double s_dot, double s_ddot) const;

protected:

  struct ViaPoint
  {
    ViaPoint(double s_, const arma::vec &pos_, double err_tol_):
    s(s_), pos(pos_), err_tol(err_tol_) {}

    double s;
    arma::vec pos;
    double err_tol;

    arma::vec z_dual;
  };
  std::list<ViaPoint> via_points;
  
  arma::mat W_mpc; ///< mpc optimized weights
  const GMP *gmp_ref;

  GMP::Ptr gmp_mpc;

  unsigned N;
  arma::rowvec dt_;

  arma::mat Aineq_slack;
  arma::mat Q_slack;
  
  arma::mat Qi;
  arma::mat QN;
  
  arma::vec Z0;
  arma::vec Z0_dual_ineq;
  arma::vec Z0_dual_eq;
  
  arma::vec pos_lb;
  arma::vec pos_ub;
  
  arma::vec vel_lb;
  arma::vec vel_ub;
  
  arma::vec accel_lb;
  arma::vec accel_ub;
  
  unsigned n_slack;
  
  bool pos_slack;
  bool vel_slack;
  bool accel_slack;
  
  arma::vec pos_slack_lim;
  arma::vec vel_slack_lim;
  arma::vec accel_slack_lim;
  
  arma::mat Phi0;
  arma::vec x0;
  
  double s_f;
  arma::mat Phi_f;
  arma::vec x_f;
  arma::vec err_tol_f;

private:

  unsigned n_dof;
  unsigned n_dof3;

  double inf;

  arma::mat I_ndof;
  arma::vec inf_ndof;
  arma::vec ones_ndof;
  arma::vec zeros_ndof;

    
}; // class GMP_MPC

} // namespace gmp_

} // namespace as64_

#endif // GMP_LIB_GMP_MPC_H