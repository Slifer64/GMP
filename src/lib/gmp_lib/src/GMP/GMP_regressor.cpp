#include <gmp_lib/GMP/GMP_regressor.h>

#include <cfloat>

namespace as64_
{

namespace gmp_
{

// public:

  GMP_regressor::GMP_regressor(unsigned N_kernels, double kernel_std_scaling)
  {
    this->kernel_fun_ptr = std::bind(&GMP_regressor::kernelFun,this,std::placeholders::_1);
    this->setKernels(N_kernels, kernel_std_scaling);
  }

  void GMP_regressor::setTruncatedKernels(bool set, double zero_tol)
  {
    // find range of x outside which psi = 0 due to finite numerical precision
    double r_min = zero_tol;
    this->x_min = this->c(0) - std::sqrt( -std::log(r_min) / this->h(0) );
    this->x_max = this->c.back() + std::sqrt( -std::log(r_min) / this->h.back() );
    
    if (set) this->kernel_fun_ptr = std::bind(&GMP_regressor::truncGaussKernel,this, std::placeholders::_1, zero_tol);
    else this->kernel_fun_ptr = std::bind(&GMP_regressor::kernelFun,this,std::placeholders::_1);
  }

  // ============================================================

  arma::vec GMP_regressor::regressVec(double x) const
  {
    arma::vec psi(this->c.size());

    // take appropriate actions when x causes phi = 0 due to finite
    // numerical precision.
    if (x < this->x_min)
    {
      psi = arma::vec().zeros(psi.size());
      psi(0) = 1;
    }
    else if (x > this->x_max)
    {
      psi = arma::vec().zeros(psi.size());
      psi.back() = 1;
    }
    else
    {
      psi = this->kernel_fun_ptr(x);
    }
    
    return psi / arma::sum(psi);
  }

  arma::vec GMP_regressor::regressVecDot(double x, double x_dot) const
  {
    if (x < this->x_min || x > this->x_max) return arma::vec().zeros(this->c.size());

    arma::vec psi = this->kernel_fun_ptr(x);
    arma::vec psi_dot = this->kernelFunDot(x, x_dot);
    double sum_psi = arma::sum(psi);
    double sum_psi_dot = arma::sum(psi_dot);

    arma::vec phi = psi / sum_psi;
    return ( psi_dot - phi*sum_psi_dot ) / sum_psi;
  }

  arma::vec GMP_regressor::regressVecDDot(double x, double x_dot, double x_ddot) const
  {
    if (x < this->x_min || x > this->x_max) return arma::vec().zeros(this->c.size());

    arma::vec psi = this->kernel_fun_ptr(x);
    arma::vec psi_dot = this->kernelFunDot(x, x_dot);
    arma::vec psi_ddot = this->kernelFunDDot(x, x_dot, x_ddot);
    double sum_psi = arma::sum(psi);
    double sum_psi_dot = arma::sum(psi_dot);
    double sum_psi_ddot = arma::sum(psi_ddot);

    arma::vec phi = psi / sum_psi;
    arma::vec phi_dot = ( psi_dot - phi*sum_psi_dot ) / sum_psi;
    return (psi_ddot - 2*phi_dot*sum_psi_dot - phi*sum_psi_ddot) / sum_psi;
  }

  arma::vec GMP_regressor::regressVec3Dot(double x, double dx, double ddx, double d3x) const
  {
    if (x < this->x_min || x > this->x_max) return arma::vec().zeros(this->c.size());

    arma::vec psi = this->kernel_fun_ptr(x);
    arma::vec psi_dot = this->kernelFunDot(x, dx);
    arma::vec psi_ddot = this->kernelFunDDot(x, dx, ddx);
    arma::vec psi_3dot = this->kernelFun3Dot(x, dx, ddx, d3x);
    double sum_psi = arma::sum(psi);
    double sum_psi_dot = arma::sum(psi_dot);
    double sum_psi_ddot = arma::sum(psi_ddot);
    double sum_psi_3dot = arma::sum(psi_3dot);

    arma::vec phi = psi / sum_psi;
    arma::vec phi_dot = ( psi_dot - phi*sum_psi_dot ) / sum_psi;
    arma::vec phi_ddot = (psi_ddot - 2*phi_dot*sum_psi_dot - phi*sum_psi_ddot) / sum_psi;
    return (psi_3dot - 3*phi_ddot*sum_psi_dot - 3*phi_dot*sum_psi_ddot - phi*sum_psi_3dot) / sum_psi;
  }

// protected:

  // =============================================================

  arma::vec GMP_regressor::kernelFun(double x) const
  {
    return arma::exp(-this->h % (arma::pow(x-this->c,2)));
  }

  arma::vec GMP_regressor::truncGaussKernel(double x, double zero_tol)
  {
    arma::vec psi = this->kernelFun(x);
    arma::uvec ind = arma::find(psi < zero_tol);
    psi.elem(ind) = arma::vec().zeros(ind.size());
    return psi;
  }

  arma::vec GMP_regressor::kernelFunDot(double x, double dx) const
  {
    arma::vec psi = this->kernel_fun_ptr(x);
    arma::vec a = (x-this->c)*dx;
    return -2*this->h%( psi % a );
  }

  arma::vec GMP_regressor::kernelFunDDot(double x, double dx, double ddx) const
  {
    arma::vec psi = this->kernel_fun_ptr(x);
    arma::vec psi_dot = this->kernelFunDot(x, dx);

    arma::vec a = (x-this->c)*dx;
    arma::vec a_dot = (x-this->c)*ddx + std::pow(dx,2);
    return -2*this->h%( psi_dot%a + psi%a_dot );
  }

  arma::vec GMP_regressor::kernelFun3Dot(double x, double dx, double ddx, double d3x) const
  {
    arma::vec psi = this->kernel_fun_ptr(x);
    arma::vec psi_dot = this->kernelFunDot(x, dx);
    arma::vec psi_ddot = this->kernelFunDDot(x, dx, ddx);

    arma::vec a = (x-this->c)*dx;
    arma::vec a_dot = (x-this->c)*ddx + std::pow(dx,2);
    arma::vec a_ddot = (x-this->c)*d3x + 3*dx*ddx;
    return -2*this->h%( psi_ddot%a + 2*psi_dot%a_dot + psi%a_ddot );
  }

  void GMP_regressor::setKernels(unsigned N_kernels, double kernel_std_scaling)
  {
    if (N_kernels < 2) throw std::runtime_error("[GMP_regressor::setKernels]: At least 2 kernels are required!");

    this->c = arma::linspace<arma::vec>(0,N_kernels-1, N_kernels)/(N_kernels-1);

    this->kernel_std_scaling = kernel_std_scaling;

    double hi = 1 / std::pow(kernel_std_scaling*(this->c(1)-this->c(0)),2);
    this->h = arma::vec().ones(N_kernels) * hi; 

    // find range of x outside which psi = 0 due to finite numerical precision
    double r_min = DBL_MAX;
    this->x_min = this->c(0) - std::sqrt( -std::log(r_min) / this->h(0) );
    this->x_max = this->c.back() + std::sqrt( -std::log(r_min) / this->h.back() );
  }

} // namespace gmp_

} // namespace as64_
