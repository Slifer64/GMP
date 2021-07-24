#include <gmp_lib/GMP/GMP_Opt.h>

#include <osqp_lib/quadprog.h>

namespace as64_
{

namespace gmp_
{

#define GMP_Opt_fun_ std::string("[GMP_Opt::") + __func__ + "]: "

  GMP_Opt::GMP_Opt(gmp_::GMP *gmp)
  {
      /*this->ex_flag_map = containers.Map('KeyType','double','ValueType','char');
      this->ex_flag_map(-10) = 'Empty... Call ''constrOpt'' first.';
      this->ex_flag_map(1) = 'Function converged to the solution x.';
      this->ex_flag_map(0) = 'Number of iterations exceeded options.MaxIterations.';

      this->ex_flag_map(-2) = 'Problem is infeasible. Or, for interior-point-convex, the step size was smaller than options.StepTolerance, but constraints were not satisfied.';
      this->ex_flag_map(-3) = 'Problem is unbounded';
      this->ex_flag_map(2) = 'Step size was smaller than options.StepTolerance, constraints were satisfied.';
      this->ex_flag_map(-6) = 'Nonconvex problem detected.';
      this->ex_flag_map(-8) = 'Unable to compute a step direction.';
      */
      this->exit_msg = "";

      this->gmp = gmp;

      this->tau = 8;

      this->setOptions(0.1,1,0, 1,1,1);
  }

  void GMP_Opt::setOptions(bool opt_pos, bool opt_vel, bool opt_accel, double pos_obj_w, double vel_obj_w, double accel_obj_w)
  {
    this->opt_pos = opt_pos;
    this->opt_vel = opt_vel;
    this->opt_accel = opt_accel;

    this->w_p = pos_obj_w;
    this->w_v = vel_obj_w;
    this->w_a = accel_obj_w;
  }

  bool GMP_Opt::optimize(unsigned num_points)
  {
    arma::rowvec x_data = arma::linspace<arma::rowvec>(0,1, num_points);
    this->optimize(x_data);
  }

  bool GMP_Opt::optimize(const arma::rowvec &x_data)
  {
    int n_ker = this->gmp->numOfKernels();
    int n_dof = this->gmp->numOfDoFs();

    double x_dot = 1/this->tau;
    double x_ddot = 0;

    bool success = true;
    this->exit_msg = "";

    // calculate cost function: J = 0.5w'Hw + f'w
    int N = x_data.size();

    arma::mat H = 1e-8*arma::eye(n_ker, n_ker); // for numerical stability
    arma::mat f = arma::mat().zeros(n_ker, n_dof);
    arma::vec phi;

    if (this->opt_pos)
    {
        arma::mat y_offset = - this->gmp->getY0() + this->gmp->getScaling()*this->gmp->getY0d();
        arma::mat H1 = arma::mat().zeros(n_ker, n_ker);
        arma::mat f1 = arma::mat().zeros(n_ker, n_dof);
        for (int i=0; i<N; i++)
        {
          phi = this->gmp->regressVec(x_data(i));
          H1 = H1 + phi*phi.t();
          f1 = f1 - phi*(this->gmp->getYd(x_data(i)) + y_offset).t();
        }
        H = H + this->w_p*H1;
        f = f + this->w_p*f1;
    }

    if (this->opt_vel)
    {
        arma::mat H2 = arma::mat().zeros(n_ker, n_ker);
        arma::mat f2 = arma::mat().zeros(n_ker, n_dof);
        for (int i=0; i<N; i++)
        {
          phi = this->gmp->regressVecDot(x_data(i), x_dot);
          H2 = H2 + phi*phi.t();
          f2 = f2 - phi*this->gmp->getYdDot(x_data(i), x_dot).t();
        }
        H = H + this->w_v*H2;
        f = f + this->w_v*f2;
    }

    if (this->opt_accel)
    {
        arma::mat H3 = arma::mat().zeros(n_ker, n_ker);
        arma::mat f3 = arma::mat().zeros(n_ker, n_dof);
        for (int i=0; i<N; i++)
        {
          phi = this->gmp->regressVecDDot(x_data(i), x_dot, x_ddot);
          H3 = H3 + phi*phi.t();
          f3 = f3 - phi*this->gmp->getYdDDot(x_data(i), x_dot, x_ddot).t();
        }
        H = H + this->w_a*H3;
        f = f + this->w_a*f3;
    }

    // inequality constraints

    arma::mat A;
    arma::mat lb;
    arma::mat ub;

    if (!this->A_p.empty())
    {
      A = arma::join_vert(A, this->A_p);
      lb = arma::join_vert(lb, this->pos_lb);
      ub = arma::join_vert(ub, this->pos_ub);
    }

    if (!this->A_v.empty())
    {
      A = arma::join_vert(A, this->A_v);
      lb = arma::join_vert(lb, this->vel_lb);
      ub = arma::join_vert(ub, this->vel_ub);
    }

    if (!this->A_a.empty())
    {
      A = arma::join_vert(A, this->A_a);
      lb = arma::join_vert(lb, this->accel_lb);
      ub = arma::join_vert(ub, this->accel_ub);
    }

    // equality constraints

    arma::mat Aeq;
    arma::mat beq;

    if (!this->Aeq_p.empty())
    {
      Aeq = arma::join_vert(Aeq, this->Aeq_p);
      beq = arma::join_vert(beq, this->pos_eq);
    }

    if (!this->Aeq_v.empty())
    {
      Aeq = arma::join_vert(Aeq, this->Aeq_v);
      beq = arma::join_vert(beq, this->vel_eq);
    }

    if (!this->Aeq_a.empty())
    {
      Aeq = arma::join_vert(Aeq, this->Aeq_a);
      beq = arma::join_vert(beq, this->accel_eq);
    }

    // H = sparse(H);
    // A = sparse(A);

    // solve optimization problem
    osqp_::QuadProgSolution solution = osqp_::quadprog(H,f, A,lb,ub, Aeq,beq);

    if (solution.success) this->gmp->W = this->gmp->getInvScaling()*solution.x.t();
    else this->exit_msg = solution.exit_msg;

    return solution.success;
  }

  void GMP_Opt::setMotionDuration(double tau)
  {
    this->tau = tau;
  }

  void GMP_Opt::setPosConstr(const arma::rowvec &x, const arma::mat &lb, const arma::mat &ub, const arma::rowvec &x_eq, const arma::mat &p_eq)
  {
    int n_ker = this->gmp->numOfKernels();
    int n_dof = this->gmp->numOfDoFs();

    arma::vec y_offset = - this->gmp->getY0() + this->gmp->getScaling()*this->gmp->getY0d();

    // Process inequality constraints
    if (!x.empty())
    {
      int m = x.size();
      if (lb.n_cols != m) throw std::runtime_error(GMP_Opt_fun_ + "Lower bounds must have " + std::to_string(m) + " columns (constraints).");
      if (ub.n_cols != m) throw std::runtime_error(GMP_Opt_fun_ + "Upper bounds must have " + std::to_string(m) + " columns (constraints).");
      if (lb.n_rows != ub.n_rows) throw std::runtime_error(GMP_Opt_fun_ + "Lower and Upper bounds must have the same number of rows (DoFs).");

      this->A_p.resize(m, n_ker);
      this->pos_lb.resize(m, n_dof);
      this->pos_ub.resize(m, n_dof);
      for (int i=0; i<m; i++)
      {
        this->pos_lb.row(i) = ( lb.col(i) + y_offset ).t();
        this->pos_ub.row(i) = ( ub.col(i) + y_offset ).t();
        this->A_p.row(i) = this->gmp->regressVec(x(i)).t();
      }
    }

    // Process equality constraints
    if (!x_eq.empty())
    {
      int m = x_eq.size();
      if (p_eq.n_cols != m) throw std::runtime_error(GMP_Opt_fun_ + "Equality values matrix must have " + std::to_string(m) + " columns (constraints).");

      this->Aeq_p.resize(m, n_ker);
      this->pos_eq.resize(m, n_dof);
      for (int i=0; i<m; i++)
      {
        this->pos_eq.row(i) = ( p_eq.col(i) + y_offset ).t();
        this->Aeq_p.row(i) = this->gmp->regressVec(x_eq(i)).t();
      }
    }

  }

  void GMP_Opt::setVelConstr(const arma::rowvec &x, const arma::mat &lb, const arma::mat &ub, const arma::rowvec &x_eq, const arma::mat &v_eq)
  {
    int n_ker = this->gmp->numOfKernels();
    double x_dot = 1 / this->tau;

    // Process inequality constraints
    if (!x.empty())
    {
      int m = x.size();
      if (lb.n_cols != m) throw std::runtime_error(GMP_Opt_fun_ + "Lower bounds must have " + std::to_string(m) + " columns (constraints).");
      if (ub.n_cols != m) throw std::runtime_error(GMP_Opt_fun_ + "Upper bounds must have " + std::to_string(m) + " columns (constraints).");
      if (lb.n_rows != ub.n_rows) throw std::runtime_error(GMP_Opt_fun_ + "Lower and Upper bounds must have the same number of rows (DoFs).");

      this->vel_lb = lb.t();
      this->vel_ub = ub.t();
      this->A_v.resize(m, n_ker);
      for (int i=0; i<m; i++) this->A_v.row(i) = this->gmp->regressVecDot(x(i), x_dot).t();

    }

    // Process equality constraints
    if (!x_eq.empty())
    {
      int m = x_eq.size();
      if (v_eq.n_cols != m) throw std::runtime_error(GMP_Opt_fun_ + "Equality values matrix must have " + std::to_string(m) + " columns (constraints).");

      this->vel_eq = v_eq.t();
      this->Aeq_v.resize(m, n_ker);
      for (int i=0; i<m; i++) this->Aeq_v.row(i) = this->gmp->regressVecDot(x_eq(i),x_dot).t();
    }

  }

  void GMP_Opt::setAccelConstr(const arma::rowvec &x, const arma::mat &lb, const arma::mat &ub, const arma::rowvec &x_eq, const arma::mat &a_eq)
  {
    int n_ker = this->gmp->numOfKernels();
    double x_dot = 1 / this->tau;
    double x_ddot = 0;

    // Process inequality constraints
    if (!x.empty())
    {
      int m = x.size();
      if (lb.n_cols != m) throw std::runtime_error(GMP_Opt_fun_ + "Lower bounds must have " + std::to_string(m) + " columns (constraints).");
      if (ub.n_cols != m) throw std::runtime_error(GMP_Opt_fun_ + "Upper bounds must have " + std::to_string(m) + " columns (constraints).");
      if (lb.n_rows != ub.n_rows) throw std::runtime_error(GMP_Opt_fun_ + "Lower and Upper bounds must have the same number of rows (DoFs).");

      this->accel_lb = lb.t();
      this->accel_ub = ub.t();
      this->A_a.resize(m, n_ker);
      for (int i=0; i<m; i++) this->A_a.row(i) = this->gmp->regressVecDDot(x(i), x_dot, x_ddot).t();
    }

    // Process equality constraints
    if (!x_eq.empty())
    {
      int m = x_eq.size();
      if (a_eq.n_cols != m) throw std::runtime_error(GMP_Opt_fun_ + "Equality values matrix must have " + std::to_string(m) + " columns (constraints).");

      this->accel_eq = a_eq.t();
      this->Aeq_a.resize(m, n_ker);
      for (int i=0; i<m; i++) this->Aeq_a.row(i) = this->gmp->regressVecDDot(x_eq(i), x_dot, x_ddot).t();

    }
  }

  void GMP_Opt::clearPosConstr()
  {
      this->A_p.clear();
      this->pos_lb.clear();
      this->pos_ub.clear();

      this->Aeq_p.clear();
      this->pos_eq.clear();
  }

  void GMP_Opt::clearVelConstr()
  {
      this->A_v.clear();
      this->vel_lb.clear();
      this->vel_ub.clear();

      this->Aeq_v.clear();
      this->vel_eq.clear();

  }

  void GMP_Opt::clearAccelConstr()
  {
      this->A_a.clear();
      this->accel_lb.clear();
      this->accel_ub.clear();

      this->Aeq_a.clear();
      this->accel_eq.clear();

  }

  void GMP_Opt::clearConstr()
  {
    this->clearPosConstr();
    this->clearVelConstr();
    this->clearAccelConstr();
  }

  std::string GMP_Opt::getExitMsg() const
  {
    return this->exit_msg;
  }

} // namespace gmp_

} // namespace as64_
