#ifndef GMP_LIB_GMP_REGRESSOR_H
#define GMP_LIB_GMP_REGRESSOR_H

#include <cmath>
#include <vector>
#include <cstring>
#include <memory>
#include <exception>
#include <fstream>
#include <armadillo>

namespace as64_
{

namespace gmp_
{

class GMP_regressor
{

public:

  /** Weighted Sum of Gaussians constructor.
  * @param[in] N_kernels: The number of kernels.
  * @param[in] kernel_std_scaling: Scaling of the kernel's std. (optional, default=1.0)
  */
  GMP_regressor(unsigned N_kernels, double kernel_std_scaling=1.0);

  // ============================================================

  /** Returns the scaled regressor vector ks*phi.
  * @param[in] x: The phase variable (must be in [0 1]).
  * @return regressor vector.
  */
  arma::vec regressVec(double x) const;

  /** Returns the scaled regressor vector 1st time derivative ks*phi_dot.
  * @param[in] x: The phase variable (must be in [0 1]).
  * @param[in] x_dot: The phase variable 1st time derivative.
  * @return regressor vector 1st time derivative.
  */
  arma::vec regressVecDot(double x, double x_dot) const;

  /** Returns the scaled regressor vector 2nd time derivative ks*phi_ddot.
  * @param[in] x: The phase variable (must be in [0 1]).
  * @param[in] x_dot: The phase variable 1st time derivative.
  * @param[in] x_ddot: The phase variable 2nd time derivative.
  * @return regressor vector 2nd time derivative.
  */
  arma::vec regressVecDDot(double x, double x_dot, double x_ddot) const;

  /** Returns the scaled regressor vector 3rd time derivative ks*phi_3dot.
  * @param[in] x: The phase variable (must be in [0 1]).
  * @param[in] x_dot: The phase variable 1st time derivative.
  * @param[in] x_ddot: The phase variable 2nd time derivative.
  * @param[in] x_3dot: The phase variable 3rd time derivative.
  * @return regressor vector 3rd time derivative.
  */
  arma::vec regressVec3Dot(double x, double dx, double ddx, double d3x) const;

protected:

  // =============================================================

  /** Returns a column vector with the values of the kernel functions.
   * @param[in] x: The phase variable.
   * @return: Column vector with the values of the kernel functions.
   */
   arma::vec kernelFun(double x) const;

   arma::vec kernelFunDot(double x, double dx) const;

   arma::vec kernelFunDDot(double x, double dx, double ddx) const;

   arma::vec kernelFun3Dot(double x, double dx, double ddx, double d3x) const;


public:

  arma::vec c; ///< N_kernels x 1 vector with the kernels' centers
  arma::vec h; ///< N_kernels x 1 vector with the kernels' inverse width

private:

  static long double zero_tol;  ///< small value used to avoid divisions with very small numbers

};

} // namespace gmp_

} // namespace as64_

#endif // GMP_LIB_GMP_REGRESSOR_H
