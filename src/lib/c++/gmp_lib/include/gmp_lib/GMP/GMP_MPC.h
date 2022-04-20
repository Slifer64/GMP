#ifndef GMP_LIB_GMP_MPC_H
#define GMP_LIB_GMP_MPC_H

// N-DoF GMP-MPC Optimization class
//

#include <gmp_lib/GMP/GMP.h>
#include <gmp_lib/math/linear_algebra.h>
#include <list>
#include <array>
#include <functional>

#include <armadillo>

#define GMP_MPC_L0G_STATISTICS

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
  
  GMP_MPC(const GMP *gmp, unsigned N_horizon, double pred_time_step, unsigned N_kernels, 
          double kernel_std_scaling, const std::array<double,3> &slack_gains, double trunc_kern_thres=1e-6);

  void setInitialState(const arma::vec &y0, const arma::vec &y0_dot, const arma::vec &y0_ddot, double s, double s_dot, double s_ddot);
  
  void setFinalState(const arma::vec &yf, const arma::vec &yf_dot, const arma::vec &yf_ddot, double s, double s_dot, double s_ddot);
  void setFinalState(const arma::vec &yf, const arma::vec &yf_dot, const arma::vec &yf_ddot, double s, double s_dot, double s_ddot, const arma::vec &err_tol);

  void setObjCostGains(double pos_gain, double vel_gain);

  void setPosLimits(const arma::vec &lb, const arma::vec &ub);
  void setPosLimits(const Eigen::VectorXd &lb, const Eigen::VectorXd &ub);
  
  void setVelLimits(const arma::vec &lb, const arma::vec &ub);
  void setVelLimits(const Eigen::VectorXd &lb, const Eigen::VectorXd &ub);
  
  void setAccelLimits(const arma::vec &lb, const arma::vec &ub);
  void setAccelLimits(const Eigen::VectorXd &lb, const Eigen::VectorXd &ub);
  
  void addViaPoint(double s, const arma::vec &y, double err_tol = 1e-6);

  void setPosSlackLimit(double s_lim);
  
  void setVelSlackLimit(double s_lim);
  
  void setAccelSlackLimit(double s_lim);

  void setCanonicalSystemFunction(std::function<std::array<double,2>(double, double)> can_sys_fun);

  GMP_MPC::Solution solve(double s, double s_dot);

  arma::vec getYd(double s) const;

  arma::vec getYdDot(double s, double s_dot) const;

  arma::vec getYdDDot(double s, double s_dot, double s_ddot) const;

  std::vector<double> elaps_t_data;
  std::vector<double> calcH_elaps_t_data;
  std::vector<double> solve_elaps_t_data;

  #ifdef GMP_MPC_L0G_STATISTICS
    Eigen::SpMat H_;
    Eigen::SpMat A_;
    double s_rec = 0.5;
  #endif

protected:

  void integratePhase(double s, double s_dot, std::vector<double> &si_data, std::vector<double> &si_dot_data, std::vector<double> &si_ddot_data);

  Eigen::SpMat getAjMat(double s, double s_dot, double s_ddot) const;

  #ifdef GMP_MPC_L0G_STATISTICS
    arma::wall_clock timer;
    arma::wall_clock timer_H;
    arma::wall_clock timer_solve;
  #endif

  std::function<std::array<double,2>(double, double)> can_sys_fun;

  struct ViaPoint
  {
    ViaPoint(double s_, const Eigen::VectorXd &pos_, double err_tol_):
    s(s_), pos(pos_), err_tol(err_tol_) {}

    double s;
    Eigen::VectorXd pos;
    Eigen::SpMat Phi;
    double err_tol;

    Eigen::VectorXd z_dual;
  };
  std::list<ViaPoint> via_points;
  int n_via; // for how many via-points to allocate space
  
  arma::mat W_mpc; ///< mpc optimized weights
  const GMP *gmp_ref;

  GMP::Ptr gmp_mpc;

  unsigned N_kernels;

  unsigned N;
  std::vector<double> dt_;

  unsigned n_ineq; // number of inequality constraints
  unsigned n_var; // number of optimization variables
  unsigned n_eq; // number of equality constraints
  
  Eigen::VectorXd Qs_vec;
  std::vector<Eigen::Triplet<double>> Aineq_slack_triplets;
  
  double sqrt_pos_gain;
  double sqrt_vel_gain;
  Eigen::VectorXd sqrt_pos_vel_gain_vec;

  // Eigen::SpMat Qi;
  // Eigen::SpMat QN;

  arma::vec Z0;
  arma::vec Z0_dual_ineq;
  arma::vec Z0_dual_eq;
  
  Eigen::VectorXd pos_lb;
  Eigen::VectorXd pos_ub;
  
  Eigen::VectorXd vel_lb;
  Eigen::VectorXd vel_ub;
  
  Eigen::VectorXd accel_lb;
  Eigen::VectorXd accel_ub;
  
  unsigned n_slack;
  
  bool pos_slack;
  bool vel_slack;
  bool accel_slack;
  
  Eigen::VectorXd pos_slack_lim;
  Eigen::VectorXd vel_slack_lim;
  Eigen::VectorXd accel_slack_lim;
  
  Eigen::SpMat Phi0;
  Eigen::VectorXd x0;
  
  double s_f;
  std::vector<Eigen::Triplet<double>> Phi_f_triplets;
  Eigen::SpMat Phi_f;
  Eigen::SpMat Qf;
  arma::vec phi_f, phi_f_dot, phi_f_ddot;
  Eigen::VectorXd x_f;
  Eigen::VectorXd err_tol_f;

private:

  unsigned n_dof;
  unsigned n_dof3;

  double inf;

  Eigen::SpMat I_ndof;
  Eigen::VectorXd inf_ndof;
  Eigen::VectorXd ones_ndof;
  Eigen::VectorXd zeros_ndof;

    
}; // class GMP_MPC

} // namespace gmp_

} // namespace as64_

#endif // GMP_LIB_GMP_MPC_H
