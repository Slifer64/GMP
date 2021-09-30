#include <osqp_lib/qp_solver.h>
#include <osqp_lib/csc_mat.h>


namespace osqp_
{

QPsolver::QPsolver(const arma::mat &H, const arma::vec &q,
         const arma::mat &A, const arma::vec &lb, const arma::vec &ub,
         const arma::mat &Aeq, const arma::vec &beq)
{
  work = NULL;
  settings = NULL;
  data = NULL;
  warm_start_x = false;
  warm_start_y = false;

  this->n_ineq = A.n_rows;
  this->n_eq = Aeq.n_rows;

  osqp_::CSC_mat P_cs(H, true);
  osqp_::CSC_mat A_cs(arma::join_vert(A, Aeq));

  this->q_ = q;
  this->lb_ = arma::join_vert(lb, beq);
  this->ub_ = arma::join_vert(ub, beq);

  settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));
  if (!settings) throw std::runtime_error("[QPsolver::QPsolver]: Failed to allocate space for OSQPSettings...");

  data = (OSQPData *)c_malloc(sizeof(OSQPData));
  if (!data) throw std::runtime_error("[QPsolver::QPsolver]: Failed to allocate space for QSQPData...");

  int n_var = A.n_cols + Aeq.n_cols;
  int n_constr = A.n_rows;

  data->n = n_var;
  data->m = n_constr;
  data->P = csc_matrix(data->n, data->n, P_cs.nnz, P_cs.dataPtr(), P_cs.rowIndPtr(), P_cs.csPtr());
  data->A = csc_matrix(data->m, data->n, A_cs.nnz, A_cs.dataPtr(), A_cs.rowIndPtr(), A_cs.csPtr());
  data->q = this->q_.memptr();
  data->l = this->lb_.memptr();
  data->u = this->ub_.memptr();

  osqp_set_default_settings(settings);
  // settings->alpha = 1.0; // Change alpha parameter
  settings->warm_start = false;
  settings->polish = 0;
  //settings->time_limit = 0;
  //settings->max_iter = 4000;
  settings->verbose = false;
  settings->warm_start = false;
}

QPsolver::~QPsolver()
{
  if (data)
  {
    if (data->A) c_free(data->A);
    if (data->P) c_free(data->P);
    c_free(data);
  }
  if (settings) c_free(settings);

  if (work) osqp_cleanup(work);
}

void QPsolver::setPrimalSolutionGuess(const arma::vec &x0)
{
  warm_start_x = true;
  this->x0 = x0;
}

void QPsolver::setDualSolutionGuess(const arma::vec &y0_ineq, const arma::vec &y0_eq)
{
  warm_start_y = true;
  this->y0_ineq = y0_ineq;
  this->y0_eq = y0_eq;
}

int QPsolver::solve()
{
  c_int exit_flag = osqp_setup(&work, data, settings);
  if (exit_flag != 0)
  {
    exit_msg = work->info->status;
    return -1;
  }

  arma::vec y0 = arma::join_vert(y0_ineq, y0_eq);
  if (warm_start_x && warm_start_y) exit_flag = osqp_warm_start(work, x0.memptr(), y0.memptr());
  else if (warm_start_x) exit_flag = osqp_warm_start_x(work, x0.memptr());
  else if (warm_start_y) exit_flag = osqp_warm_start_y(work, y0.memptr());

  if (exit_flag != 0)
  {
    exit_msg = work->info->status;
    return -1;
  }

  // Solve Problem
  osqp_solve(work);

  exit_msg = work->info->status;

  c_int sol_status = work->info->status_val;

  if (sol_status == -3 || sol_status == -4 || sol_status == -7 || sol_status == -10) return -1;

  x0 = arma::mat(work->solution->x, data->n, 1, true);
  y0_ineq = arma::mat(work->solution->y, n_ineq, 1, true);
  y0_eq = arma::mat(work->solution->y+n_ineq, n_eq, 1, true);

  if (sol_status == 1) return 0;
  else return 1;
}

arma::vec QPsolver::getPrimalSolution() const
{
  return x0;
}

arma::vec QPsolver::getIneqDualSolution() const
{
  return y0_ineq;
}

arma::vec QPsolver::getEqDualSolution() const
{
  return y0_eq;
}

std::string QPsolver::getExitMsg() const
{
  return exit_msg;
}

} // namespace osqp_

