#ifndef GMP_LIB_GMP_NDOF_UPDATE_H
#define GMP_LIB_GMP_NDOF_UPDATE_H

// N-DoF GMP Update class
//

#include <gmp_lib/GMP/GMP_nDoF.h>

namespace as64_
{

namespace gmp_
{

class GMP_nDoF_Update
{

public:

  typedef std::shared_ptr<GMP_nDoF_Update> Ptr;

  /**GMP constructor.
   *  @param[in] gmp: n_DoF dmp.
   */
  GMP_nDoF_Update(gmp_::GMP_nDoF *gmp);

  void initSigmaw();

  void initSigmaWfromMsr(const arma::rowvec &x_data);

  void enableSigmawUpdate(bool flag);

  void setSigmaW(const arma::mat &Sw);

  arma::mat getSigmaW() const;

  double setMsrNoiseVar(double rv);

  // ==================================================
  // ===============   Online update  =================
  // ==================================================

  void updatePos(double x, const arma::vec &y, double r_n=-1);

  void updateVel(double x, double x_dot, const arma::vec &y_dot, double r_n=-1);

  void updateAccel(double x, double x_dot, double x_ddot, const arma::vec &y_ddot, double r_n=-1);


  /** Updates the weights based on 'n' measurements.
   *  @param[in] s: 1 x n vector of type gmp_::Phase, where the j-th entry is the phase for the j-th measurement.
   *  @param[in] Z: n_dof x n matrix, where the j-th column is the j-th measurement.
   *  @param[in] type: 1 x n vector of type gmp_::UPDATE_TYPE, where the j-th entry has the update type for the j-th measurement.
   *  @param[in] r_n: 1 x n vector, where the j-th entry is the noise variance for the j-th measurement (optional, default=this->rv).
   */
  void updateWeights(const std::vector<gmp_::Phase> &s, arma::mat Z, const std::vector<gmp_::UPDATE_TYPE> &type, arma::rowvec r_n={});

private:

  gmp_::GMP_nDoF *gmp; ///< pointer to nDoF GMP object

  arma::mat Sigma_w; ///< weights covariance for trajectory update
  double rv; ///< noise variance for trajectory update
  bool enable_Sigma_w_update;

}; // GMP_nDoF_Update

} // namespace gmp_

} // namespace as64_

#endif // GMP_LIB_GMP_NDOF_UPDATE_H
