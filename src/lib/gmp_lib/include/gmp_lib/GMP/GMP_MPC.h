#ifndef GMP_LIB_GMP_MPC_H
#define GMP_LIB_GMP_MPC_H

// N-DoF GMP-MPC Optimization class
//

#include <gmp_lib/GMP/GMP.h>


#include <cfloat>
#include <armadillo>

#include <osqp_lib/constants.h>


namespace as64_
{

namespace gmp_
{

class GMP_MPC
{

public:
  
  GMP_MPC(const GMP *gmp, unsigned N_horizon, double pred_time_step, unsigned N_kernels, double kernel_std_scaling, const std::vector<bool> &slack_flags, const std::vector<double> &slack_gains)
  {
    this->INF = OSQP_INFTY;

    unsigned n_dof = gmp.numOfDoFs();
    
    this->gmp_ref = gmp;
    
    this->gmp_mpc.reset(new GMP(n_dof, N_kernels, kernel_std_scaling));
    this->gmp_mpc->setTruncatedKernels(true,1e-8);
    
    this->N = N_horizon;
    this->dt_ = pred_time_step*arma::rowvec().ones(this->N);

    bool opt_pos = 1;
    bool opt_vel = 0;
    this->pos_slack = slack_flags[0];
    this->vel_slack = slack_flags[1];
    this->accel_slack = slack_flags[2];
    this->n_slack = this->pos_slack + this->vel_slack + this->accel_slack;

    // this->Aineq_slack = [];
    // this->Q_slack = [];
    if (this->pos_slack)
    {
      this->Q_slack = blkdiag(this->Q_slack, slack_gains[0]); 
      this->Aineq_slack = arma::join_horiz(this->Aineq_slack, arma::join_vert(-arma::vec().ones(n_dof), arma::vec().zeros(2*n_dof)) );
    }
    if (this->vel_slack)
    {
      this->Q_slack = blkdiag(this->Q_slack, slack_gains[1]);
      this->Aineq_slack = arma::join_horiz(this->Aineq_slack, arma::join_vert(arma::vec().zeros(n_dof), -arma::vec().ones(n_dof), arma::vec().zeros(n_dof) ) );
    }
    if (this->accel_slack)
    {
      this->Q_slack = blkdiag(this->Q_slack, slack_gains[2]);
      this->Aineq_slack = arma::join_horiz(this->Aineq_slack, arma::join_vert(arma::vec().zeros(2*n_dof), -arma::vec().ones(n_dof)) );
    }
    //this->Q_slack = sparse(this->Q_slack);
    //this->Aineq_slack = sparse(this->Aineq_slack);

    // State tracking gains: (x(i) - xd(i)).t()*Qi*(x(i) - xd(i))
    arma::mat I_ndof = arma::mat().eye(n_dof,n_dof);
    this->Qi = blkdiag( opt_pos*I_ndof , opt_vel*10*I_ndof);
    this->QN = this->Qi; // blkdiag( 100*speye(n_dof,n_dof) , 1*speye(n_dof,n_dof));

    unsigned n_dof3 = 3*n_dof;
    this->Z0 = arma::vec().zeros(n_dof*N_kernels + this->n_slack);
    this->Z0_dual_ineq = aram::vec().zeros(this->N*n_dof3 + this->n_slack);
    this->Z0_dual_eq = arma::vec().zeros(n_dof3);
    
    arma::vec inf_ndof = this->INF*arma::vec().ones(n_dof);
    this->setPosLimits(-inf_ndof, inf_ndof));
    this->setVelLimits(-inf_ndof, inf_ndof);
    this->setAccelLimits(-inf_ndof, inf_ndof);
    
    this->setPosSlackLimit(this->INF);
    this->setVelSlackLimit(this->INF);
    this->setAccelSlackLimit(this->INF);
  }

  arma::mat blkdiag(const arma::mat &A, const arma::mat &B)
  {
    arma::mat C(A.n_rows + B.n_rows, A.n_cols + B.n_cols);
    C.submat(0,0, A.n_rows-1, A.n_cols-1) = A;
    C.submat(A.n_rows, A.n_cols, C.n_rows-1, C.n_cols-1) = B;
  }

  arma::mat blkdiag(const arma::mat &A, double b)
  {
    arma::mat C(A.n_rows + 1, A.n_cols + 1);
    C.submat(0,0, A.n_rows-1, A.n_cols-1) = A;
    C.submat(A.n_rows, A.n_cols, C.n_rows-1, C.n_cols-1) = b;
  }

  
  void setInitialState(const arma::vec &y0, const arma::vec &y0_dot, const arma::vec &y0_ddot, double s, double s_dot, double s_ddot)
  {
    unsigned n_dof = this->gmp_mpc->numOfDoFs();
    arma::vec phi0 = this->gmp_mpc->regressVec(s);
    arma::vec phi0_dot = this->gmp_mpc->regressVecDot(s, s_dot);
    arma::vec phi0_ddot = this->gmp_mpc->regressVecDDot(s, s_dot, s_ddot);

    arma::mat In = arma::mat().eye(n_dof, n_dof);
    this->Phi0 = arma::join_vert( arma::kron(In, phi0.t()), arma::kron(In, phi0_dot.t()), arma::kron(In, phi0_ddot.t()) );
    this->x0 = arma::join_vert(y0, y0_dot, y0_ddot);
  }
  
  void setFinalState(const arma::vec &yf, const arma::vec &yf_dot, const arma::vec &yf_ddot, double s, double s_dot, double s_ddot)
  {
    unsigned n_dof = this->gmp_mpc->numOfDoFs();
    this->s_f = s;
    arma::vec phi_f = this->gmp_mpc->regressVec(s);
    arma::vec phi_f_dot = this->gmp_mpc->regressVecDot(s, s_dot);
    arma::vec phi_f_ddot = this->gmp_mpc->regressVecDDot(s, s_dot, s_ddot);

    arma::mat In = arma::mat().eye(n_dof, n_dof;
    this->Phi_f = arma::join_vert( arma::kron(In, phi_f.t()), arma::kron(In, phi_f_dot.t()), arma::kron(In, phi_f_ddot.t()) );
    this->x_f = arma::join_vert(yf, yf_dot, yf_ddot);
  }

  void setPosLimits(const arma::vec &lb, const arma::vec &ub)
  {
    this->pos_lb = lb;
    this->pos_ub = ub;
  }
  
  void setVelLimits(const arma::vec &lb, const arma::vec &ub)
  {
    this->vel_lb = lb;
    this->vel_ub = ub;
  }
  
  void setAccelLimits(const arma::vec &lb, const arma::vec &ub)
  {   
    this->accel_lb = lb;
    this->accel_ub = ub;
  }
  
  void setPosSlackLimit(double s_lim)
  {   
    this->pos_slack_lim = s_lim;
  }
  
  void setVelSlackLimit(double s_lim)
  {    
    this->vel_slack_lim = s_lim;
  }
  
  void setAccelSlackLimit(double s_lim)
  {   
    this->accel_slack_lim = s_lim;
  }

  void [exit_flag, y, y_dot, y_ddot, slack_var] = solve(s, s_dot, s_ddot)
  {  
      y = [];
      y_dot = []; 
      y_ddot = [];
      slack_var = [];
      
      N_kernels = this->gmp_mpc->numOfKernels();
      n_dof = this->gmp_ref.numOfDoFs();
      n_dof3 = 3*n_dof; % for pos, vel, accel
      
      n = n_dof * N_kernels + this->n_slack;
      
      H = sparse(n,n);
      q = zeros(n,1);
      % add slacks to bounds. I.e. one could optionaly define with slacks
      % bounds the maximum allowable error
      Aineq = sparse(this->N*n_dof3 + this->n_slack, n);
      Aineq(end-this->n_slack+1:end,end-this->n_slack+1:end) = speye(this->n_slack);

      % DMP phase variable
      si = s;
      si_dot = s_dot;
      si_ddot = s_ddot;

      for i=1:this->N

          yd_i = this->gmp_ref.getYd(si);
          dyd_i = this->gmp_ref.getYdDot(si, si_dot);

          phi = this->gmp_mpc->regressVec(si);
          phi_dot = this->gmp_mpc->regressVecDot(si, si_dot);
          phi_ddot = this->gmp_mpc->regressVecDDot(si, si_dot, si_ddot);

          if (i==this->N), Qi_ = this->QN;
          else, Qi_ = this->Qi;

          Psi = [kron(speye(n_dof),phi.t()); kron(speye(n_dof),phi_dot.t())];
          xd_i = [yd_i; dyd_i];

          H = H + blkdiag(Psi.t()*Qi_*Psi, this->Q_slack);
          q = q - [Psi.t()*Qi_*xd_i; zeros(this->n_slack,1)];

          Aineq_i = [kron(speye(n_dof),phi.t()); kron(speye(n_dof),phi_dot.t()); kron(speye(n_dof),phi_ddot.t())];
          Aineq((i-1)*n_dof3+1 : i*n_dof3, :) = [Aineq_i, this->Aineq_slack];

          si = si + si_dot*this->dt_(i);
          si_dot = si_dot + si_ddot*this->dt_(i);
          % si_ddot = ... (if it changes too)
      }

      H = 1e-6*speye(n) + (H+H.t())/2; % to account for numerical errors

      % Here we set the slacks to be unbounded, but in general, we could
      % specify as bounds the maximum allowable violation from the
      % kinematic bounds.
      z_min = [this->pos_lb; this->vel_lb; this->accel_lb];
      z_max = [this->pos_ub; this->vel_ub; this->accel_ub];
      
      slack_lim = [];
      if (this->pos_slack), slack_lim = [slack_lim; this->pos_slack_lim];
      if (this->vel_slack), slack_lim = [slack_lim; this->vel_slack_lim];
      if (this->accel_slack), slack_lim = [slack_lim; this->accel_slack_lim];
      
      Z_min = [repmat(z_min, this->N,1); -slack_lim];
      Z_max = [repmat(z_max, this->N,1); slack_lim];

      Aeq = [this->Phi0, sparse(n_dof3, this->n_slack)];
      beq = [this->x0];

      if (si >= this->s_f)
          Aeq = [Aeq; [this->Phi_f, sparse(n_dof3, this->n_slack)] ];
          beq = [beq; this->x_f];

          n_eq_plus = size(Aeq,1) - length(this->Z0_dual_eq);
          if (n_eq_plus>0), this->Z0_dual_eq = [this->Z0_dual_eq; zeros(n_eq_plus, 1)];
      }

      %% ===========  solve optimization problem  ==========

      A_osqp = [Aineq; Aeq];
      lb = [Z_min; beq];
      ub = [Z_max; beq];

      Z0_dual = [this->Z0_dual_ineq; this->Z0_dual_eq];

      % Create an OSQP object
      osqp_solver = osqp;
      //osqp_solver.setup(H, q, A_osqp, lb, ub, 'warm_start',false, 'verbose',false, 'eps_abs',1e-4, 'eps_rel',1e-5);%, 'max_iter',20000);
      osqp_solver.warm_start('x', this->Z0, 'y',Z0_dual);

      res = osqp_solver.solve();

      this->exit_msg = ""; % clear exit message
      exit_flag = 0;
      if ( res.info.status_val ~= 1)
          %res.info
          this->exit_msg = res.info.status;
          exit_flag = 1;
          if (res.info.status_val == -3 || res.info.status_val == -4 || res.info.status_val == -7 || res.info.status_val == -10)
              exit_flag = -1;
              return; 
          }
      }

      Z = res.x;
      this->Z0 = Z;
      Z0_dual = res.y;
      n_ineq = size(Aineq,1);
      this->Z0_dual_ineq = Z0_dual(1:n_ineq);
      this->Z0_dual_eq = Z0_dual(n_ineq+1:end);
      
      w = Z(1:end-this->n_slack);
      slack_var = Z(end-this->n_slack+1:end);
      W = reshape(w, N_kernels, n_dof).t();

      %% --------  Generate output  --------

      y = W*this->gmp_mpc->regressVec(s);
      y_dot = W*this->gmp_mpc->regressVecDot(s, s_dot);
      y_ddot = W*this->gmp_mpc->regressVecDDot(s, s_dot, s_ddot);

      this->setInitialState(y, y_dot, y_ddot, s, s_dot, s_ddot);
      
  }
  
  std::string getExitMsg() const
  {  
    return this->exit_msg;
      
  {


protected:

  std::string exit_msg;
  
  const GMP *gmp_ref;

  GMP::Ptr gmp_mpc;

  unsigned N;
  arma::rowvec<double> dt_;

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
  
  double pos_slack_lim;
  double vel_slack_lim;
  double accel_slack_lim;
  
  arma::mat Phi0;
  arma::vec x0;
  
  double s_f;
  arma::mat Phi_f;
  arma::vec x_f;

private:

  double INF;

    
}; // class GMP_MPC

} // namespace gmp_

} // namespace as64_

#endif // GMP_LIB_GMP_MPC_H
