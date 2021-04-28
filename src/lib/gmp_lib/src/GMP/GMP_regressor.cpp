#include <gmp_lib/GMP/GMP_regressor.h>

namespace as64_
{

namespace gmp_
{

long double GMP_regressor::zero_tol = 1e-200;

// public:

  GMP_regressor::GMP_regressor(unsigned N_kernels, double kernel_std_scaling)
  {
    this->c = arma::linspace<arma::vec>(0,N_kernels-1, N_kernels)/(N_kernels-1);
    double hi = 1 / std::pow(kernel_std_scaling*(this->c(1)-this->c(0)),2);
    this->h = arma::vec().ones(N_kernels) * hi;
  }

  // ============================================================

  arma::vec GMP_regressor::regressVec(double x) const
  {
    arma::vec psi = this->kernelFun(x);
    return psi / (arma::sum(psi) + this->zero_tol);
  }

  arma::vec GMP_regressor::regressVecDot(double x, double x_dot) const
  {
    arma::vec psi = this->kernelFun(x);
    arma::vec psi_dot = this->kernelFunDot(x, x_dot);
    double sum_psi = arma::sum(psi);
    double sum_psi_dot = arma::sum(psi_dot);

    arma::vec phi = psi / ( arma::sum(sum_psi) + this->zero_tol );
    return ( psi_dot - phi*sum_psi_dot ) / ( sum_psi + this->zero_tol);
  }

  arma::vec GMP_regressor::regressVecDDot(double x, double x_dot, double x_ddot) const
  {
    arma::vec psi = this->kernelFun(x);
    arma::vec psi_dot = this->kernelFunDot(x, x_dot);
    arma::vec psi_ddot = this->kernelFunDDot(x, x_dot, x_ddot);
    double sum_psi = arma::sum(psi);
    double sum_psi_dot = arma::sum(psi_dot);
    double sum_psi_ddot = arma::sum(psi_ddot);

    arma::vec phi = psi / ( arma::sum(sum_psi) + this->zero_tol );
    arma::vec phi_dot = ( psi_dot - phi*sum_psi_dot ) / ( sum_psi + this->zero_tol);
    return (psi_ddot - 2*phi_dot*sum_psi_dot - phi*sum_psi_ddot) / ( sum_psi + this->zero_tol);
  }

  arma::vec GMP_regressor::regressVec3Dot(double x, double dx, double ddx, double d3x) const
  {
    arma::vec psi = this->kernelFun(x);
    arma::vec psi_dot = this->kernelFunDot(x, dx);
    arma::vec psi_ddot = this->kernelFunDDot(x, dx, ddx);
    arma::vec psi_3dot = this->kernelFun3Dot(x, dx, ddx, d3x);
    double sum_psi = arma::sum(psi);
    double sum_psi_dot = arma::sum(psi_dot);
    double sum_psi_ddot = arma::sum(psi_ddot);
    double sum_psi_3dot = arma::sum(psi_3dot);

    arma::vec phi = psi / ( arma::sum(sum_psi) + this->zero_tol );
    arma::vec phi_dot = ( psi_dot - phi*sum_psi_dot ) / ( sum_psi + this->zero_tol);
    arma::vec phi_ddot = (psi_ddot - 2*phi_dot*sum_psi_dot - phi*sum_psi_ddot) / ( sum_psi + this->zero_tol);
    return (psi_3dot - 3*phi_ddot*sum_psi_dot - 3*phi_dot*sum_psi_ddot - phi*sum_psi_3dot) / ( sum_psi + this->zero_tol);
  }

// protected:

  // =============================================================

   arma::vec GMP_regressor::kernelFun(double x) const
   {
     return arma::exp(-this->h % (arma::pow(x-this->c,2)));
   }

   arma::vec GMP_regressor::kernelFunDot(double x, double dx) const
   {
     arma::vec psi = this->kernelFun(x);
     arma::vec a = (x-this->c)*dx;
     return -2*this->h%( psi % a );
   }

   arma::vec GMP_regressor::kernelFunDDot(double x, double dx, double ddx) const
   {
     arma::vec psi = this->kernelFun(x);
     arma::vec psi_dot = this->kernelFunDot(x, dx);

     arma::vec a = (x-this->c)*dx;
     arma::vec a_dot = (x-this->c)*ddx + std::pow(dx,2);
     return -2*this->h%( psi_dot%a + psi%a_dot );
   }

   arma::vec GMP_regressor::kernelFun3Dot(double x, double dx, double ddx, double d3x) const
   {
     arma::vec psi = this->kernelFun(x);
     arma::vec psi_dot = this->kernelFunDot(x, dx);
     arma::vec psi_ddot = this->kernelFunDDot(x, dx, ddx);

     arma::vec a = (x-this->c)*dx;
     arma::vec a_dot = (x-this->c)*ddx + std::pow(dx,2);
     arma::vec a_ddot = (x-this->c)*d3x + 3*dx*ddx;
     return -2*this->h%( psi_ddot%a + 2*psi_dot%a_dot + psi%a_ddot );
   }

} // namespace gmp_

} // namespace as64_
