#ifndef OSQP_QP_SOLVER_H
# define OSQP_QP_SOLVER_H

#include <memory>
#include <osqp_lib/osqp.h>
#include <osqp_lib/csc_mat.h>

#include <armadillo>

#ifndef NULL_MAT
  #define NULL_MAT arma::mat(0,0)
#endif

#ifndef NULL_VEC
  #define NULL_VEC arma::vec(0)
#endif

namespace osqp_
{

class QPsolver
{
public:
  QPsolver(const arma::mat &H, const arma::vec &q,
           const arma::mat &A, const arma::vec &lb, const arma::vec &ub,
           const arma::mat &Aeq, const arma::vec &beq);

  ~QPsolver();

  void setPrimalSolutionGuess(const arma::vec &x0);

  void setDualSolutionGuess(const arma::vec &y0_ineq, const arma::vec &y0_eq);

  int solve();

  arma::vec getPrimalSolution() const;

  arma::vec getIneqDualSolution() const;

  arma::vec getEqDualSolution() const;

  std::string getExitMsg() const;

  OSQPSettings  *settings;

private:

  std::shared_ptr<osqp_::CSC_mat> P_cs;
  std::shared_ptr<osqp_::CSC_mat> A_cs;

  OSQPWorkspace *work;
  OSQPData      *data;

  int n_ineq;
  int n_eq;

  arma::vec q_;
  arma::vec lb_;
  arma::vec ub_;

  bool warm_start_x;
  arma::vec x0;

  bool warm_start_y;
  arma::vec y0_ineq, y0_eq;

  std::string exit_msg;

}; // class QPsolver


} // namespace osqp_

#endif // OSQP_QP_SOLVER_H
