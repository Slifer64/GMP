#ifndef GMP_LIB_GMP_NDOF_H
#define GMP_LIB_GMP_NDOF_H

// N-DoF GMP class
// Generalized movement primitive.
//

#include <gmp_lib/GMP/GMP_regressor.h>
#include <gmp_lib/utils.h>
#include <gmp_lib/TrajScale/TrajScale_Prop.h>
#include <gmp_lib/TrajScale/TrajScale_Rot_min.h>
#include <gmp_lib/TrajScale/TrajScale_Rot_wb.h>

namespace as64_
{

namespace gmp_
{

class GMP_nDoF_Update; // forward declaration
class GMP_nDoF_IO; // forward declaration

class GMP_nDoF: public GMP_regressor
{

public:

  typedef std::shared_ptr<GMP_nDoF> Ptr;

  /** GMP constructor.
   *  @param[in] n_dofs: number of degrees of freedom.
   *  @param[in] N_kernels: the number of kernels
   *  @param[in] kern_std_scale: Scaling for std of kernels (optional, default=1).
   */
  GMP_nDoF(unsigned n_dofs, unsigned N_kernels, double kern_std_scale=1.0);


  /** Returns the number of DoFs.
   *  return: number of DoFs.
   */
  unsigned numOfDoFs() const;

  /** Returns the number of kernels.
   *  @return: number of kernels.
   */
  unsigned numOfKernels() const;

  /** Sets the initial position.
   *  @param[in] y0: initial position.
   */
  void setY0(const arma::vec &y0);

  arma::vec getY0() const;

  arma::vec getY0d();


  /** Set goal position.
   *  @param[in] g: goal position.
   */
  void setGoal(const arma::vec &g);

  arma::vec getGoal() const;

  arma::mat getScaling() const;

  arma::mat getInvScaling() const;

  /** Set the trajectory spatial scaling method.
   *  @param[in] traj_scale_obj: pointer to an object of type @TrajScale.
   */
  void setScaleMethod(const std::shared_ptr<gmp_::TrajScale> &traj_scale_obj);

  /** Trains the GMP.
   * @param[in] train_method: the training method to use, as a string ('LWR', 'LS').
   * @param[in] x: Row vector with the canonical timestamps (in [0 1]) of the training data points.
   * @param[in] yd_data: Matrix with the desired potition for each DoF in each row.
   * @param[out] train_error: The training error expressed as the mse error.
   * \note The timestamps-data need not be sequencial temporarily.
   */
  void train(const std::string &train_method, const arma::rowvec &x, const arma::mat &yd_data, arma::vec *train_err=NULL, arma::mat *Sw=NULL);

  /** Calculates the time derivatives of the GMP's states.
   * @param[in] s: Object of type @GMP_phase.
   * @param[in] y: 'y' state of the GMP.
   * @param[in] z: 'z' state of the GMP.
   * @param[in] y_c: coupling term for the dynamical equation of the 'y' state (optional, default=0).
   * @param[in] z_c: coupling term for the dynamical equation of the 'z' state (optional, default=0).
   */
  void update(const gmp_::Phase &s, const arma::vec &y, const arma::vec &z, arma::vec y_c={0}, arma::vec z_c={0});

  /** Returns the 'y' state time derivative.
   * Call after @update.
   * @return: time derivative of 'y' state.
   */
  arma::vec getYdot() const;


  /** Returns the 'z' state time derivative.
   * Call after @update.
   * @return: time derivative of 'z' state.
   */
  arma::vec getZdot() const;

  /** Returns the GMP's acceleration.
   * Call after @update.
   * @param[in] yc_dot: time derivative of 'y' state coupling (optional, default=0).
   * @return: acceleration.
   */
  arma::vec getYddot(arma::vec yc_dot={0}) const;


  /** Calculates the GMP's acceleration.
   * @param[in] s: Object of type @GMP_phase.
   * @param[in] y: 'y' state of the GMP.
   * @param[in] y_dot: time derivative of 'y' state.
   * @param[in] y_c: coupling term for the dynamical equation of the 'y' state (optional, default=0).
   * @param[in] z_c: coupling term for the dynamical equation of the 'z' state (optional, default=0).
   * @param[in] yc_dot: time derivative of the 'y' state coupling (optional, default=0).
   * @return: acceleration.
   */
  arma::vec calcYddot(const gmp_::Phase &s, const arma::vec &y, const arma::vec &y_dot, arma::vec y_c={0}, arma::vec z_c={0}, arma::vec yc_dot={0});

  /** Returns the scaled desired position.
   * @param[in] x: phase variable.
   * @return: scaled desired position.
   */
  arma::vec getYd(double x) const;


  /** Returns the scaled desired velocity.
   * @param[in] x: phase variable.
   * @param[in] x_dot: 1st time derivative of the phase variable.
   * @return: scaled desired velocity.
   */
  arma::vec getYdDot(double x, double x_dot) const;


  /** Returns the scaled desired acceleration.
   * @param[in] x: phase variable.
   * @param[in] x_dot: 1st time derivative of the phase variable.
   * @param[in] x_ddot: 2nd time derivative of the phase variable.
   * @return: scaled desired acceleration.
   */
  arma::vec getYdDDot(double x, double x_dot, double x_ddot=0) const;


public: // properties

  gmp_::TrajScale::Ptr traj_sc; ///< object of type @TrajScale

  // weights
  arma::mat W; ///< num_DoF x num_Kernels matrix where each row contrains the weights for each DoF

  // impedance params
  arma::vec K; ///< num_DoF x 1 stiffness vector
  arma::vec D; ///< num_DoF x 1 stiffness vector

protected: // properties

  friend class GMP_nDoF_Update;
  friend class GMP_nDoF_IO;

  arma::vec Y0; ///< initial position
  arma::vec Yg; ///< target position

  arma::vec Y0d; ///< trained initial position
  arma::vec Ygd; ///< trained target position

  // output state
  arma::vec y_dot; ///< position derivative
  arma::vec z_dot; ///< scaled velocity derivative

}; // GMP_nDoF

} // namespace gmp_

} // namespace as64_

#endif // GMP_LIB_GMP_NDOF_H
