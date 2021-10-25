#include <gmp_lib/GMP/GMP_MPC.h>

#include <cfloat> // for infinity

//#include <osqp_lib/osqp.h>
//#include <osqp_lib/constants.h>
//#include <osqp_lib/csc_mat.h>
#include <osqp_lib/qp_solver.h>

namespace as64_
{

namespace gmp_
{


arma::mat blkdiag(const arma::mat &A, const arma::mat &B)
{
  arma::mat C(A.n_rows + B.n_rows, A.n_cols + B.n_cols);
  C.submat(0,0, A.n_rows-1, A.n_cols-1) = A;
  C.submat(A.n_rows, A.n_cols, C.n_rows-1, C.n_cols-1) = B;
  return C;
}

arma::mat blkdiag(const arma::mat &A, double b)
{
  arma::mat C(A.n_rows + 1, A.n_cols + 1);
  if (!A.empty()) C.submat(0,0, A.n_rows-1, A.n_cols-1) = A;
  C.submat(A.n_rows, A.n_cols, C.n_rows-1, C.n_cols-1) = b;
  return C;
}


GMP_MPC::GMP_MPC(const GMP *gmp, unsigned N_horizon, double pred_time_step, unsigned N_kernels, double kernel_std_scaling, const std::vector<double> &slack_gains)
{
  this->inf = OSQP_INFTY;

  this->settings.time_limit = 0;
  this->settings.max_iter = 4000;
  this->settings.abs_tol = 1e-4;
  this->settings.rel_tol = 1e-5;

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

  this->pos_slack = slack_gains[0] > 0;
  this->vel_slack = slack_gains[1] > 0;
  this->accel_slack = slack_gains[2] > 0;
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

  // State tracking gains: (x(i) - xd(i)).t()*Qi*(x(i) - xd(i))
  setObjCostGains(1.0, 0.0);
  // this->Qi = blkdiag( opt_pos_gain*I_ndof, opt_vel_gain*I_ndof );
  // this->QN = this->Qi; //blkdiag( 100*I_ndof, 1*I_ndof );

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

void GMP_MPC::setObjCostGains(double pos_gain, double vel_gain)
{
  this->Qi = blkdiag( pos_gain*I_ndof, vel_gain*I_ndof );
  this->QN = this->Qi; //blkdiag( 100*I_ndof, 1*I_ndof );
}

void GMP_MPC::setInitialState(const arma::vec &y0, const arma::vec &y0_dot, const arma::vec &y0_ddot, double s, double s_dot, double s_ddot)
{
  arma::vec phi0 = this->gmp_mpc->regressVec(s);
  arma::vec phi0_dot = this->gmp_mpc->regressVecDot(s, s_dot);
  arma::vec phi0_ddot = this->gmp_mpc->regressVecDDot(s, s_dot, s_ddot);

  this->Phi0 = arma::join_vert( arma::kron(I_ndof, phi0.t()), arma::kron(I_ndof, phi0_dot.t()), arma::kron(I_ndof, phi0_ddot.t()) );
  this->x0 = arma::join_vert(y0, y0_dot, y0_ddot);
}

void GMP_MPC::setFinalState(const arma::vec &yf, const arma::vec &yf_dot, const arma::vec &yf_ddot, double s, double s_dot, double s_ddot)
{
  setFinalState(yf, yf_dot, yf_ddot, s, s_dot, s_ddot, arma::vec().zeros(yf.size()) );
}

void GMP_MPC::setFinalState(const arma::vec &yf, const arma::vec &yf_dot, const arma::vec &yf_ddot, double s, double s_dot, double s_ddot, const arma::vec &err_tol)
{
  this->s_f = s;
  this->err_tol_f = arma::kron( err_tol, arma::vec().ones(n_dof,1) );
  arma::vec phi_f = this->gmp_mpc->regressVec(s);
  arma::vec phi_f_dot = this->gmp_mpc->regressVecDot(s, s_dot);
  arma::vec phi_f_ddot = this->gmp_mpc->regressVecDDot(s, s_dot, s_ddot);

  this->Phi_f = arma::join_vert( arma::kron(I_ndof, phi_f.t()), arma::kron(I_ndof, phi_f_dot.t()), arma::kron(I_ndof, phi_f_ddot.t()) );
  this->x_f = arma::join_vert(yf, yf_dot, yf_ddot);
}

void GMP_MPC::setPosLimits(const arma::vec &lb, const arma::vec &ub)
{
  this->pos_lb = lb;
  this->pos_ub = ub;
}

void GMP_MPC::setVelLimits(const arma::vec &lb, const arma::vec &ub)
{
  this->vel_lb = lb;
  this->vel_ub = ub;
}

void GMP_MPC::setAccelLimits(const arma::vec &lb, const arma::vec &ub)
{   
  this->accel_lb = lb;
  this->accel_ub = ub;
}

void GMP_MPC::setPosSlackLimit(double s_lim)
{   
  this->pos_slack_lim = s_lim;
}

void GMP_MPC::setVelSlackLimit(double s_lim)
{    
  this->vel_slack_lim = s_lim;
}

void GMP_MPC::setAccelSlackLimit(double s_lim)
{   
  this->accel_slack_lim = s_lim;
}

GMP_MPC::Solution GMP_MPC::solve(double s, double s_dot, double s_ddot)
{
  arma::vec y = x0.subvec(0, n_dof-1);
  arma::vec y_dot = x0.subvec(n_dof, 2*n_dof-1);
  arma::vec y_ddot = x0.subvec(2*n_dof, 3*n_dof-1);
  arma::vec slack_var = arma::vec().zeros(n_slack);

  unsigned N_kernels = this->gmp_mpc->numOfKernels();
  unsigned n = n_dof * N_kernels + this->n_slack;

  arma::mat H = 1e-6*arma::mat().eye(n,n); // for numerical stability
  arma::vec q = arma::vec().zeros(n);
  // add slacks to bounds. I.e. one could optionaly define with slacks
  // bounds the maximum allowable error
  arma::mat Aineq = arma::mat().zeros(this->N*n_dof3 + this->n_slack, n);
  if (this->n_slack)
  {
    int i_end = Aineq.n_rows - 1;
    int j_end = Aineq.n_cols - 1;
    Aineq.submat(i_end-this->n_slack+1, j_end-this->n_slack+1, i_end, j_end) = arma::mat().eye(this->n_slack,this->n_slack);
  }

  arma::vec lb_i = arma::join_vert(this->pos_lb, this->vel_lb, this->accel_lb);
  arma::vec ub_i = arma::join_vert(this->pos_ub, this->vel_ub, this->accel_ub);

  arma::vec slack_lim;
  if (this->pos_slack) slack_lim = arma::join_vert(slack_lim, arma::vec({this->pos_slack_lim}) );
  if (this->vel_slack) slack_lim = arma::join_vert(slack_lim, arma::vec({this->vel_slack_lim}) );
  if (this->accel_slack) slack_lim = arma::join_vert(slack_lim, arma::vec({this->accel_slack_lim}) );

  arma::vec lb = arma::join_vert( arma::repmat(lb_i, this->N,1), -slack_lim);
  arma::vec ub = arma::join_vert( arma::repmat(ub_i, this->N,1), slack_lim);

  // DMP phase variable
  double si = s;
  double si_dot = s_dot;
  double si_ddot = s_ddot;

  arma::mat Hi(n,n);
  arma::vec qi(n);
  if (this->n_slack) Hi.submat(n-n_slack, n-n_slack, n-1, n-1) = this->Q_slack;

  for (int i=1; i<=this->N; i++)
  {
    arma::vec yd_i = this->gmp_ref->getYd(si);
    arma::vec dyd_i = this->gmp_ref->getYdDot(si, si_dot);

    if (si >= 1)
    {
      yd_i = x_f.subvec(0,n_dof-1);
      dyd_i = x_f.subvec(n_dof, 2*n_dof-1);
    }

    arma::vec phi = this->gmp_mpc->regressVec(si);
    arma::vec phi_dot = this->gmp_mpc->regressVecDot(si, si_dot);
    arma::vec phi_ddot = this->gmp_mpc->regressVecDDot(si, si_dot, si_ddot);

    arma::mat Qi_;
    if (i==this->N) Qi_ = this->QN;
    else Qi_ = this->Qi;

    arma::mat Psi = arma::join_vert( arma::kron(I_ndof,phi.t()), arma::kron(I_ndof,phi_dot.t()) );
    arma::vec xd_i = arma::join_vert(yd_i, dyd_i);

    Hi.submat(0,0, n-n_slack-1, n-n_slack-1) = Psi.t()*Qi_*Psi;
    H = H + Hi;
    qi.subvec(0, n-n_slack-1) = -Psi.t()*Qi_*xd_i;
    q = q + qi;

    arma::mat Aineq_i = arma::join_vert( arma::kron(I_ndof,phi.t()), arma::kron(I_ndof,phi_dot.t()), arma::kron(I_ndof,phi_ddot.t()) );
    Aineq.rows((i-1)*n_dof3, i*n_dof3-1) = arma::join_horiz(Aineq_i, this->Aineq_slack);

    si = si + si_dot*this->dt_(i-1);
    si_dot = si_dot + si_ddot*this->dt_(i-1);
    // si_ddot = ... (if it changes too)
  }

  arma::mat Aeq = arma::join_horiz( this->Phi0, arma::mat().zeros(n_dof3, this->n_slack) );
  arma::vec beq = this->x0;

  // if (si >= this->s_f)
  // {
    Aineq = arma::join_vert(Aineq, arma::join_horiz(this->Phi_f, arma::mat().zeros(n_dof3, this->n_slack) ) );
    lb = arma::join_vert(lb, this->x_f - this->err_tol_f);
    ub = arma::join_vert(ub, this->x_f + this->err_tol_f);

    int n_ineq_plus = Aineq.n_rows - this->Z0_dual_ineq.size();
    if (n_ineq_plus>0) this->Z0_dual_ineq = arma::join_vert(this->Z0_dual_ineq, arma::vec().zeros(n_ineq_plus) );

    // Aeq = arma::join_vert(Aeq, arma::join_horiz(this->Phi_f, arma::mat().zeros(n_dof3, this->n_slack) ) );
    // beq = arma::join_vert(beq, this->x_f);

    // int n_eq_plus = Aeq.n_rows - this->Z0_dual_eq.size();
    // if (n_eq_plus>0) this->Z0_dual_eq = arma::join_vert(this->Z0_dual_eq, arma::vec().zeros(n_eq_plus) );
  // }

  // ===========  solve optimization problem  ==========
  osqp_::QPsolver solver(H, q, Aineq, lb, ub, Aeq, beq);
  solver.setPrimalSolutionGuess(this->Z0);
  solver.setDualSolutionGuess(this->Z0_dual_ineq, this->Z0_dual_eq);

  solver.settings->time_limit = this->settings.time_limit;
  solver.settings->max_iter = this->settings.max_iter;
  solver.settings->eps_abs = this->settings.abs_tol;
  solver.settings->eps_rel = this->settings.rel_tol;
  solver.settings->warm_start = false;
  solver.settings->verbose = false;

  int exit_flag = solver.solve();
  std::string exit_msg = solver.getExitMsg();
  if (exit_flag >= 0)
  {
    arma::vec X = solver.getPrimalSolution();
    arma::vec w = X.subvec(0, X.n_elem-1-this->n_slack);
    if (this->n_slack) slack_var = X.subvec(X.n_elem-this->n_slack, X.n_elem-1);
    arma::mat W = arma::reshape(w, N_kernels, n_dof).t();

    this->Z0 = X;
    this->Z0_dual_ineq = solver.getIneqDualSolution();
    this->Z0_dual_eq = solver.getEqDualSolution();

    // --------  Generate output  --------
    y = W*this->gmp_mpc->regressVec(s);
    y_dot = W*this->gmp_mpc->regressVecDot(s, s_dot);
    y_ddot = W*this->gmp_mpc->regressVecDDot(s, s_dot, s_ddot);

    if (this->n_slack)
    {
      // check just in case...
      if ( !arma::find( arma::abs(slack_var) > slack_lim ).is_empty() )
      {
        exit_msg = "Slack variable limits violated";
        exit_flag = -1;
        goto end_solve;  
      }
    }

    this->setInitialState(y, y_dot, y_ddot, s, s_dot, s_ddot);
  }

  end_solve:
  Solution sol;
  sol.y = y;
  sol.y_dot = y_dot;
  sol.y_ddot = y_ddot;
  sol.slack_var = slack_var;
  sol.exit_flag = exit_flag;
  sol.exit_msg = exit_msg;
  return sol;
}
  
} // namespace gmp_

} // namespace as64_
