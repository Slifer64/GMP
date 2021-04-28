#ifndef AS64_GMP_ORIENT_H_
#define AS64_GMP_ORIENT_H_

// GMPo class
// For encoding Cartesian Orientation.
//

#include <gmp_lib/GMP/GMP.h>
#include <gmp_lib/GMP/GMP_nDoF.h>
#include <gmp_lib/utils.h>

namespace as64_
{

namespace gmp_
{


class GMPo : public GMP_nDoF
{
// ===================================
// =======  Public Functions  ========
// ===================================
public:

  /* Constructs a GMP defined in the quat-log space.
   * @param[in] N_kernels: vector with the number of kernels for each dim of the quat log.
   * @param[in] D: vector with the damping for each dim of the quat log.
   * @param[in] K: vector with the stiffness for each dim of the quat log.
   * @param[in] kernels_std_scaling: Scaling for std of kernels (optional, default=2).
   * \note: Each of the arguments 'N_kernels', 'D', 'K' can be 1x1 or a 3x1 vector.
   *        If a 1x1 vector is passed, the same value is passed to all three dim of quat log.
   */
  GMPo(arma::uvec N_kernels, arma::vec D, arma::vec K, double kernels_std_scaling=1.0);


  /* Trains the GMPo.
   * @param[in] train_method: the training method to use, as a string ('LWR', 'LS').
   * @param[in] Time: Row vector with the timestamps of the training data points.
   * @param[in] Quat_data: Matrix with the desired unit quaternion (4x1 vector) in each column.
   * @param[out] train_error: The training error expressed as the mse error.
   */
  void train(const std::string &train_method, const arma::rowvec &Time, const arma::mat &Quat_data, arma::vec *train_error=0);


  // See @GMP_nDoF::update
  // void update(const gmp_::Phase &s, const arma::vec &y, const arma::vec &z, arma::vec y_c=arma::vec({0}), arma::vec z_c=arma::vec({0}));


  // See @GMP_nDoF.getYdot
  // arma::vec getYdot() const;


  // See @GMP_nDoF.getYdot
  // arma::vec getZdot() const;


  // See @GMP_nDoF.getYdot
  // arma::vec getYddot(arma::vec yc_dot=arma::vec({0})) const;


  // See @GMP_nDoF.getYdot
  // arma::vec calcYddot(const gmp_::Phase &s, const arma::vec &y, const arma::vec &y_dot, arma::vec yc=arma::vec({0}), arma::vec zc=arma::vec({0}), arma::vec yc_dot=arma::vec({0})) const;


  // See @GMP_nDoF.getYdot
  // unsigned numOfKernels(int i) const;


  /* Sets the initial orientation.
   * @param[in] Q0: Initial orientation (as unit quaternion).
   */
  void setQ0(const arma::vec &Q0);

  arma::vec getQ0() const { return Q0; }

  arma::vec getQd0() const { return Qd0; }

  arma::vec getQdf() const { return q2quat(getYdf(),Qd0); }


  /* Sets goal/target orientation.
   * @param[in] Qg: Goal/target orientation (as unit quaternion).
   */
  void setQg(const arma::vec &Qg);


  // See @GMP_nDoF.deepCopy
  // cp_obj = deepCopy()

  arma::vec getQd(double x) const;
  arma::vec getRotVeld(double x, double x_dot) const;
  arma::vec getRotAcceld(double x, double x_dot, double x_ddot) const;


  // See @GMP_nDoF.getYd
  // arma::vec getYd(double x) const;


  // See @GMP_nDoF.getYdDot
  // arma::vec getYdDot(double x, double x_dot) const;


  // See @GMP_nDoF.getYdDDot
  // arma::vec getYdDDot(double x, double x_dot, double x_ddot) const;

  arma::vec getRotVel(const arma::vec &Q) const;

  arma::vec getRotAccel(const arma::vec &Q, double tau_dot=0, const arma::vec &Yc_dot=arma::vec().zeros(3)) const;

  arma::vec calcRotAccel(const gmp_::Phase &s, const arma::vec &Q, const arma::vec &rotVel, const arma::vec &Qg,
                         const arma::vec &Yc=arma::vec().zeros(3), const arma::vec &Zc=arma::vec().zeros(3), const arma::vec &Yc_dot=arma::vec().zeros(3)) const;

  arma::vec getY(const arma::vec &Q) const;

  arma::vec getZ(const arma::vec &rotVel, const arma::vec &Q) const;


  /* Export the GMP model to a file.
   * @param[in] filename: The name of the file.
   */
  void exportToFile(const std::string &filename) const;

  /* Import a GMP model from a file.
   * @param[in] filename: The name of the file.
   */
  static std::shared_ptr<GMPo> importFromFile(const std::string &filename);

  /* Write the GMP model to a file.
   * @param[in] fid: Object of type @FileIO associated with the file.
   * @param[in] prefix: The prefix that will be used for writing the names of all GMP params (optional, default="").
   */
  void writeToFile(FileIO &fid, const std::string &prefix="") const;

  /* Reads the GMP model from a file.
   * @param[in] fid: Object of type @FileIO associated with the file.
   * @param[in] prefix: The prefix that will be used for reading the names of all GMP params (optional, default="").
   */
  void readFromFile(FileIO &fid, const std::string &prefix="");


// ===========================================================
// ***********************************************************
// ************      Static Public Functions      ************
// ***********************************************************
// ===========================================================

  /* Expresses a given quaternion w.r.t. the initial orientation.
   * @param[in] Q: Orientation as unit quaternion.
   * @param[in] Q0: Initial quaternion.
   * @return: The orientation w.r.t. Q0, i.e. Q1 = Q*Q0^{-1}.
   */
  static arma::vec quatTf(const arma::vec &Q, const arma::vec &Q0);


  /* Returns the log of a given orientation w.r.t. the initial orientation.
   * @param[in] Q: Orientation as unit quaternion.
   * @param[in] Q0: Initial quaternion.
   * @return: The logarithm of the Q w.r.t. Q0, i.e. q = log(Q*Q0^{-1}).
   */
  static arma::vec quat2q(const arma::vec &Q, const arma::vec &Q0);


  /* Returns the quaternion Q given the initial orientation Q0 and the log of Q w.r.t. Q0.
   * @param[in] q: Logarithm of orientation w.r.t. the initial orientation.
   * @param[in] Q0: Initial orientation.
   * @return: The orientation corresponding to log, i.e. Q = exp(q)*Q0
   */
  static arma::vec q2quat(const arma::vec &q, const arma::vec &Q0);


  /* Returns derivative of log given the rotational velocity and orientation (expressed w.r.t. the initial orientation)
   * @param[in] rotVel: Rotational velocity.
   * @param[in] Q1: Orientation expressed w.r.t. the initial orientation.
   * @return: Derivative of log.
   */
  static arma::vec rotVel2qdot(const arma::vec &rotVel, const arma::vec &Q1);


  /** \brief TODO doc.*/
  static arma::vec qdot2rotVel(const arma::vec &qdot, const arma::vec &Q1);


  /** \brief TODO doc.*/
  static arma::vec rotAccel2qddot(const arma::vec &rotAccel, const arma::vec &rotVel, const arma::vec &Q1);


  /** \brief TODO doc.*/
  static arma::vec qddot2rotAccel(const arma::vec &ddeo, const arma::vec &rotVel, const arma::vec &Q1);


  /** \brief TODO doc.*/
  static arma::mat jacobQq(const arma::vec &Q1);


  /** \brief TODO doc.*/
  static arma::mat jacobqQ(const arma::vec &Q1);


  /** \brief TODO doc.*/
  static arma::mat jacobDotqQ(const arma::vec &Q1, const arma::vec &rotVel);


  /** \brief TODO doc.*/
  static arma::mat jacobDotQq(const arma::vec &Q1, const arma::vec &rotVel);

// =======================================
// =======  Protected Properties  ========
// =======================================
protected:

  arma::vec Q0; ///< initial orientation (as unit quaternion)

  arma::vec Qd0; ///< initial demonstrated orientation (as unit quaternion)


  // =======  Static properties  ========

  static double zero_tol;

};

} // namespace gmp_

} // namespace as64_

#endif // AS64_GMP_ORIENT_H_
