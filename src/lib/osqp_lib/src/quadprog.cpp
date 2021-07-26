#include <osqp_lib/quadprog.h>

namespace osqp_
{

template<typename T>
T *vector_deepCopy(const T *src, c_int n)
{
  T *dest = (T *)c_malloc(n*sizeof(T));
  memcpy( (void *)dest, (void *)src, sizeof(T)*n);
  return dest;
}

OSQPPolish *OSQPPolish_deepCopy(const OSQPPolish *obj)
{
  return const_cast<OSQPPolish *>(obj);
  
  OSQPPolish *cp = (OSQPPolish *)c_malloc(sizeof(OSQPPolish));

  return cp;
}

LinSysSolver *LinSysSolver_deepCopy(const LinSysSolver *obj)
{
  return const_cast<LinSysSolver *>(obj);

  LinSysSolver *cp = (LinSysSolver *)c_malloc(sizeof(LinSysSolver));

  return cp;
}

OSQPTimer *OSQPTimer_deepCopy(const OSQPTimer *obj)
{
  return const_cast<OSQPTimer *>(obj);

  OSQPTimer *cp = (OSQPTimer *)c_malloc(sizeof(OSQPTimer));

  return cp;
}



OSQPData *OSQPData_deepCopy(const OSQPData *obj)
{
  return const_cast<OSQPData *>(obj);

  OSQPData *cp = (OSQPData *)c_malloc(sizeof(OSQPData));

  return cp;
}

OSQPSettings *OSQPSettings_deepCopy(const OSQPSettings *obj)
{
  return const_cast<OSQPSettings *>(obj);

  OSQPSettings *cp = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));

  return cp;
}

OSQPScaling *OSQPScaling_deepCopy(const OSQPScaling *obj)
{
  return const_cast<OSQPScaling *>(obj);

  OSQPScaling *cp = (OSQPScaling *)c_malloc(sizeof(OSQPScaling));

  return cp;
}


OSQPSolution *OSQPSolution_deepCopy(const OSQPSolution *obj)
{
  return const_cast<OSQPSolution *>(obj);

  OSQPSolution *cp = (OSQPSolution *)c_malloc(sizeof(OSQPSolution));

  return cp;
}


OSQPInfo *OSQPInfo_deepCopy(const OSQPInfo *obj)
{
  return const_cast<OSQPInfo *>(obj);

  OSQPInfo *cp = (OSQPInfo *)c_malloc(sizeof(OSQPInfo));

  return cp;
}


OSQPWorkspace *OSQPWorkspace_deepCopy(const OSQPWorkspace *work)
{
  OSQPWorkspace *cp = (OSQPWorkspace *)c_malloc(sizeof(OSQPWorkspace));

  cp->data = OSQPData_deepCopy(work->data);

  c_int n = cp->data->n;
  c_int m = cp->data->m;

  cp->linsys_solver = LinSysSolver_deepCopy(work->linsys_solver);

  # ifndef EMBEDDED
    cp->pol = OSQPPolish_deepCopy(work->pol);
  # endif // ifndef EMBEDDED

  // Vector used to store a vectorized rho parameter
  cp->rho_vec = vector_deepCopy(work->rho_vec, m); // vector of rho values
  cp->rho_inv_vec = vector_deepCopy(work->rho_inv_vec, m); // vector of inv rho values

  # if EMBEDDED != 1
    cp->constr_type = vector_deepCopy(work->constr_type, m); // Type of constraints: loose (-1), equality (1), inequality (0)
  # endif // if EMBEDDED != 1

  // Iterates
  cp->x = vector_deepCopy(work->x, n); // Iterate x
  cp->y = vector_deepCopy(work->x, m); // Iterate y
  cp->z = vector_deepCopy(work->x, m); // Iterate z

  cp->xz_tilde = vector_deepCopy(work->xz_tilde, n+m); // Iterate xz_tilde
  
  cp->x_prev = vector_deepCopy(work->x_prev, n); // Previous x

  cp->z_prev = vector_deepCopy(work->z_prev, m); // Previous z

  // Primal and dual residuals workspace variables
  cp->Ax = vector_deepCopy(work->Ax, m); // scaled A * x
  cp->Px = vector_deepCopy(work->Px, n); // scaled P * x
  cp->Aty = vector_deepCopy(work->Aty, n); // scaled A' * y

  //Primal infeasibility variables
  cp->delta_y = vector_deepCopy(work->delta_y, m); // difference between consecutive dual iterates
  cp->Atdelta_y = vector_deepCopy(work->Atdelta_y, n); // A' * delta_y

  // Dual infeasibility variables
  cp->delta_x = vector_deepCopy(work->delta_x, n); // difference between consecutive primal iterates
  cp->Pdelta_x = vector_deepCopy(work->Pdelta_x, n); // P * delta_x
  cp->Adelta_x = vector_deepCopy(work->Adelta_x, m); // A * delta_x

  //Temporary vectors used in scaling
  cp->D_temp = vector_deepCopy(work->D_temp, n); // temporary primal variable scaling vectors
  cp->D_temp_A = vector_deepCopy(work->D_temp_A, n); // temporary primal variable scaling vectors storing norms of A columns
  cp->E_temp = vector_deepCopy(work->E_temp, m); // temporary constraints scaling vectors storing norms of A' columns


  cp->settings = OSQPSettings_deepCopy(work->settings); // problem settings
  cp->scaling = OSQPScaling_deepCopy(work->scaling); // scaling vectors
  cp->solution = OSQPSolution_deepCopy(work->solution); // problem solution
  cp->info = OSQPInfo_deepCopy(work->info); // solver information

  # ifdef PROFILING
    cp->timer = OSQPTimer_deepCopy(work->timer);
    cp->first_run = work->first_run;
    cp->clear_update_time = work->clear_update_time;
    cp->first_run = work->first_run;
    cp->rho_update_from_solve = work->rho_update_from_solve;
  # endif // ifdef PROFILING

  # ifdef PRINTING
    cp->summary_printed = work->summary_printed;
  # endif // ifdef PRINTING

  return cp;
}

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
