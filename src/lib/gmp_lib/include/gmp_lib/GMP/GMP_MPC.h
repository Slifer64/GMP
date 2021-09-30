#ifndef GMP_LIB_GMP_MPC_H
#define GMP_LIB_GMP_MPC_H

// N-DoF GMP-MPC Optimization class
//

#include <gmp_lib/GMP/GMP.h>


#include <cfloat>
#include <armadillo>

#include <osqp_lib/osqp.h>
#include <osqp_lib/constants.h>
#include <osqp_lib/csc_mat.h>

#include <ros/ros.h>

namespace as64_
{

namespace gmp_
{

class GMP_MPC
{

public:

  struct Solution
  {
    arma::vec y;
    arma::vec y_dot;
    arma::vec y_ddot;

    arma::vec slack_var;

    int exit_flag;
    std::string exit_msg;
  };
  
  GMP_MPC(const GMP *gmp, unsigned N_horizon, double pred_time_step, unsigned N_kernels, double kernel_std_scaling, const std::vector<bool> &slack_flags, const std::vector<double> &slack_gains)
  {
    this->inf = OSQP_INFTY;

    this->n_dof = gmp->numOfDoFs();
    this->n_dof3 = 3*n_dof;
    this->I_ndof = arma::mat().eye(n_dof, n_dof);
    this->inf_ndof = this->inf*arma::vec().ones(n_dof);
    this->ones_ndof = arma::vec().ones(n_dof);
    this->zeros_ndof = arma::vec().zeros(n_dof);
     
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
      this->Aineq_slack = arma::join_horiz(this->Aineq_slack, arma::join_vert( -ones_ndof, arma::vec().zeros(2*n_dof) ) );
    }
    if (this->vel_slack)
    {
      this->Q_slack = blkdiag(this->Q_slack, slack_gains[1]);
      this->Aineq_slack = arma::join_horiz(this->Aineq_slack, arma::join_vert( zeros_ndof, -ones_ndof, zeros_ndof ) );
    }
    if (this->accel_slack)
    {
      this->Q_slack = blkdiag(this->Q_slack, slack_gains[2]);
      this->Aineq_slack = arma::join_horiz(this->Aineq_slack, arma::join_vert( arma::vec().zeros(2*n_dof), -ones_ndof ) );
    }
    //this->Q_slack = sparse(this->Q_slack);
    //this->Aineq_slack = sparse(this->Aineq_slack);

    // State tracking gains: (x(i) - xd(i)).t()*Qi*(x(i) - xd(i))
    this->Qi = blkdiag( opt_pos*I_ndof, opt_vel*10*I_ndof );
    this->QN = blkdiag( 100*I_ndof, 1*I_ndof ); // this->Qi;

    this->Z0 = arma::vec().zeros(n_dof*N_kernels + this->n_slack);
    this->Z0_dual_ineq = arma::vec().zeros(this->N*n_dof3 + this->n_slack);
    this->Z0_dual_eq = arma::vec().zeros(n_dof3);
    
    this->setPosLimits(-inf_ndof, inf_ndof);
    this->setVelLimits(-inf_ndof, inf_ndof);
    this->setAccelLimits(-inf_ndof, inf_ndof);
    
    this->setPosSlackLimit(this->inf);
    this->setVelSlackLimit(this->inf);
    this->setAccelSlackLimit(this->inf);
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
    arma::vec phi0 = this->gmp_mpc->regressVec(s);
    arma::vec phi0_dot = this->gmp_mpc->regressVecDot(s, s_dot);
    arma::vec phi0_ddot = this->gmp_mpc->regressVecDDot(s, s_dot, s_ddot);

    this->Phi0 = arma::join_vert( arma::kron(I_ndof, phi0.t()), arma::kron(I_ndof, phi0_dot.t()), arma::kron(I_ndof, phi0_ddot.t()) );
    this->x0 = arma::join_vert(y0, y0_dot, y0_ddot);
  }
  
  void setFinalState(const arma::vec &yf, const arma::vec &yf_dot, const arma::vec &yf_ddot, double s, double s_dot, double s_ddot)
  {
    this->s_f = s;
    arma::vec phi_f = this->gmp_mpc->regressVec(s);
    arma::vec phi_f_dot = this->gmp_mpc->regressVecDot(s, s_dot);
    arma::vec phi_f_ddot = this->gmp_mpc->regressVecDDot(s, s_dot, s_ddot);

    this->Phi_f = arma::join_vert( arma::kron(I_ndof, phi_f.t()), arma::kron(I_ndof, phi_f_dot.t()), arma::kron(I_ndof, phi_f_ddot.t()) );
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

  GMP_MPC::Solution solve(double s, double s_dot, double s_ddot)
  {  

    arma::vec y;
    arma::vec y_dot; 
    arma::vec y_ddot;
    arma::vec slack_var;

    unsigned N_kernels = this->gmp_mpc->numOfKernels();

    unsigned n = n_dof * N_kernels + this->n_slack;

    arma::mat H = 1e-6*arma::mat().eye(n,n); // for numerical stability
    arma::vec q = arma::vec().zeros(n);
    // add slacks to bounds. I.e. one could optionaly define with slacks
    // bounds the maximum allowable error
    arma::mat Aineq = arma::mat().zeros(this->N*n_dof3 + this->n_slack, n);
    int i_end = Aineq.n_rows - 1;
    int j_end = Aineq.n_cols - 1;
    Aineq.submat(i_end-this->n_slack+1, j_end-this->n_slack+1, i_end, j_end) = arma::mat().eye(this->n_slack,this->n_slack);

    arma::vec z_min = arma::join_vert(this->pos_lb, this->vel_lb, this->accel_lb);
    arma::vec z_max = arma::join_vert(this->pos_ub, this->vel_ub, this->accel_ub);

    arma::vec slack_lim;
    if (this->pos_slack) slack_lim = arma::join_vert(slack_lim, arma::vec({this->pos_slack_lim}) );
    if (this->vel_slack) slack_lim = arma::join_vert(slack_lim, arma::vec({this->vel_slack_lim}) );
    if (this->accel_slack) slack_lim = arma::join_vert(slack_lim, arma::vec({this->accel_slack_lim}) );

    arma::vec Z_min = arma::join_vert( arma::repmat(z_min, this->N,1), -slack_lim);
    arma::vec Z_max = arma::join_vert( arma::repmat(z_max, this->N,1), slack_lim);

    // DMP phase variable
    double si = s;
    double si_dot = s_dot;
    double si_ddot = s_ddot;

    for (int i=1; i<=this->N; i++)
    {
      arma::vec yd_i = this->gmp_ref->getYd(si);
      arma::vec dyd_i = this->gmp_ref->getYdDot(si, si_dot);

      arma::vec phi = this->gmp_mpc->regressVec(si);
      arma::vec phi_dot = this->gmp_mpc->regressVecDot(si, si_dot);
      arma::vec phi_ddot = this->gmp_mpc->regressVecDDot(si, si_dot, si_ddot);

      arma::vec Qi_;
      if (i==this->N) Qi_ = this->QN;
      else Qi_ = this->Qi;

      arma::mat Psi = arma::join_vert( arma::kron(I_ndof,phi.t()), arma::kron(I_ndof,phi_dot.t()) );
      arma::vec xd_i = arma::join_vert( yd_i, dyd_i);

      H = H + blkdiag(Psi.t()*Qi_*Psi, this->Q_slack);
      q = q - arma::join_vert( Psi.t()*Qi_*xd_i, arma::vec().zeros(this->n_slack) );

      arma::mat Aineq_i = arma::join_vert( arma::kron(I_ndof,phi.t()), arma::kron(I_ndof,phi_dot.t()), arma::kron(I_ndof,phi_ddot.t()) );
      Aineq.rows((i-1)*n_dof3, i*n_dof3-1) = arma::join_horiz(Aineq_i, this->Aineq_slack);

      si = si + si_dot*this->dt_(i);
      si_dot = si_dot + si_ddot*this->dt_(i);
      // si_ddot = ... (if it changes too)
    }

    arma::mat Aeq = arma::join_horiz( this->Phi0, arma::mat().zeros(n_dof3, this->n_slack) );
    arma::vec beq = this->x0;

    if (si >= this->s_f)
    {
      Aeq = arma::join_vert(Aeq, arma::join_horiz(this->Phi_f, arma::mat().zeros(n_dof3, this->n_slack) ) );
      beq = arma::join_vert(beq, this->x_f);

      int n_eq_plus = Aeq.n_rows - this->Z0_dual_eq.size();
      if (n_eq_plus>0) this->Z0_dual_eq = arma::join_vert(this->Z0_dual_eq, arma::vec().zeros(n_eq_plus) );
    }

    // ===========  solve optimization problem  ==========

    arma::mat A_osqp = arma::join_vert(Aineq, Aeq);
    arma::vec lb = arma::join_vert(Z_min, beq);
    arma::vec ub = arma::join_vert(Z_max, beq);

    arma::vec Z0_dual = arma::join_vert(this->Z0_dual_ineq, this->Z0_dual_eq);

    std::string exit_msg = "";
    arma::vec X, Y;
    int exit_flag = step_solve(H, q, A_osqp, lb, ub, this->Z0, Z0_dual, &X, &Y, &exit_msg);

    this->Z0 = X;
    Z0_dual = Y;
    unsigned n_ineq = Aineq.n_rows;
    this->Z0_dual_ineq = Z0_dual.subvec(0, n_ineq-1);
    this->Z0_dual_eq = Z0_dual.subvec(n_ineq, Z0_dual.n_elem-1);

    arma::vec w = X.subvec(0, X.n_elem-1-this->n_slack);
    if (this->n_slack) slack_var = X(X.n_elem-this->n_slack, X.n_elem-1);
    arma::mat W = arma::reshape(w, N_kernels, n_dof).t();

    // --------  Generate output  --------

    y = W*this->gmp_mpc->regressVec(s);
    y_dot = W*this->gmp_mpc->regressVecDot(s, s_dot);
    y_ddot = W*this->gmp_mpc->regressVecDDot(s, s_dot, s_ddot);

    this->setInitialState(y, y_dot, y_ddot, s, s_dot, s_ddot);

    Solution sol;
    sol.y = y;
    sol.y_dot = y_dot;
    sol.y_ddot = y_ddot;
    sol.slack_var = slack_var;
    sol.exit_flag = exit_flag;
    sol.exit_msg = exit_msg;
    return sol;
  }


protected:

    int step_solve(const arma::mat &H, arma::mat &q, const arma::mat &A, arma::vec &lb, arma::vec &ub, const arma::vec &X0, const arma::vec &Y0,
                    arma::vec *X, arma::vec *Y, std::string *exit_msg=0)
    {
        int n_var = A.n_cols;
        int n_constr = A.n_rows;

        osqp_::CSC_mat P_cs(H, true);
        osqp_::CSC_mat A_cs(A);

        // Exitflag
        c_int exitflag = 0;

        // Workspace structures
        OSQPWorkspace *work;
        OSQPSettings  *settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));
        OSQPData      *data     = (OSQPData *)c_malloc(sizeof(OSQPData));

        // Populate data
        if (data)
        {
            data->n = n_var;
            data->m = n_constr;
            data->P = csc_matrix(data->n, data->n, P_cs.nnz, P_cs.dataPtr(), P_cs.rowIndPtr(), P_cs.csPtr());
            data->A = csc_matrix(data->m, data->n, A_cs.nnz, A_cs.dataPtr(), A_cs.rowIndPtr(), A_cs.csPtr());
            data->q = &(q(0));
            data->l = &(lb(0));
            data->u = &(ub(0));
        }
        else throw std::runtime_error("Failed to allocate space for QSQPData...");

        // Define solver settings as default
        if (settings)
        {
            osqp_set_default_settings(settings);
            // printQSQPSettings(settings);
            settings->alpha = 1.0; // Change alpha parameter
            settings->warm_start = false;
            settings->polish = 0;
            //settings->time_limit = 0;
            //settings->max_iter = 4000;
            settings->verbose = false;
        }
        else throw std::runtime_error("Failed to allocate space for OSQPSettings...");

        // ========  Setup workspace  =========
         exitflag = osqp_setup(&work, data, settings);
         if (exitflag != 0) throw std::runtime_error(std::string("[osqp_setup]: ") + work->info->status + "\n");

        exitflag = osqp_warm_start(work, X0.memptr(), Y0.memptr());
        if (exitflag != 0) throw std::runtime_error(std::string("[osqp_setup]: ") + work->info->status + "\n");

        // Solve Problem
        osqp_solve(work);

        c_int sol_status = work->info->status_val;
        int ret;

        if (sol_status != 1)
        {
            if (exit_msg) *exit_msg = work->info->status;
            ret = 1;
            if (sol_status == -3 || sol_status == -4 || sol_status == -7 || sol_status == -10) ret = -1;
        }

        *X = arma::mat(work->solution->x, n_var, 1, true);
        *Y = arma::mat(work->solution->y, n_var, 1, true);

        // Cleanup
        if (data)
        {
            if (data->A) c_free(data->A);
            if (data->P) c_free(data->P);
            c_free(data);
        }
        if (settings) c_free(settings);

        return ret;
    }
  
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
  
  double pos_slack_lim;
  double vel_slack_lim;
  double accel_slack_lim;
  
  arma::mat Phi0;
  arma::vec x0;
  
  double s_f;
  arma::mat Phi_f;
  arma::vec x_f;

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
