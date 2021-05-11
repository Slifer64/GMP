#include <osqp_lib/quadprog.h>

namespace osqp_
{

QuadProgSolution quadprog(const arma::mat &H, const arma::mat &f, const arma::mat &A, const arma::mat &lb, const arma::mat &ub,
              const arma::mat &Aeq, const arma::mat &beq)
{
  arma::mat f_ = f;
  arma::mat lb_ = arma::join_vert(lb, beq);
  arma::mat ub_ = arma::join_vert(ub, beq);

  int n_var = A.n_cols;
  int n_constr = A.n_rows + Aeq.n_rows;

  CSC_mat P_cs(H, true);
  CSC_mat A_cs(arma::join_vert(A, Aeq));

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
  }

  // Define solver settings as default
  if (settings)
  {
    osqp_set_default_settings(settings);
    settings->alpha = 1.0; // Change alpha parameter
    settings->verbose = false;
  }

  int n_dim = f.n_cols;
  arma::mat W(n_var, n_dim);

  data->q = &(f_.col(0)(0));
  data->l = &(lb_.col(0)(0));
  data->u = &(ub_.col(0)(0));

  // Setup workspace
  exitflag = osqp_setup(&work, data, settings);
  if (exitflag != 0) throw std::runtime_error(std::string("[osqp_setup]: ") + work->info->status + "\n");

  QuadProgSolution solution;
  solution.x.resize(n_var, n_dim);
  solution.success = true;
  solution.exit_msg = "";

  for (int k=0; k<n_dim; k++)
  {
//    // Setup workspace
//    exitflag = osqp_setup(&work, data, settings);
//    exitflag = osqp_warm_start_x(work, &(x0.col(k)(0)) );

    // Solve Problem
    osqp_solve(work);

    if (work->info->status_val != OSQP_SOLVED)
    {
      // std::string err_msg = work->info->status;
      // std::cerr << "Dim " << (k+1) << ": " << err_msg << "\n";
      if (!solution.exit_msg.empty()) solution.exit_msg += "\n";
      solution.exit_msg = std::string("Dim ") + std::to_string(k+1) + ": " + work->info->status;
      solution.success = false;
    }
    else solution.x.col(k) = arma::mat(work->solution->x, n_var, 1, true);

    if (k != n_dim-1)
    {
      osqp_update_lin_cost(work, &(f_.col(k + 1)(0)));
      osqp_update_bounds(work, &(lb_.col(k + 1)(0)), &(ub_.col(k + 1)(0)));
    }
  }
  solution.success = true;

  // Cleanup
  if (data)
  {
    if (data->A) c_free(data->A);
    if (data->P) c_free(data->P);
    c_free(data);
  }
  if (settings) c_free(settings);

  return solution;
}

} // namespace osqp_
