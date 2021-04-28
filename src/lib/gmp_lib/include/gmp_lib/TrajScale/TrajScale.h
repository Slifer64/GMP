#ifndef GMP_TRAJ_SCALE_H
#define GMP_TRAJ_SCALE_H

// Trajectory scaling class
// Calculates the scaling matrix for scaling spatially a trajectory from a
// nominal start/final position {Y0d, Ygd} to a new start/final position {Y0, Yg}.
// The scaling is performed according to the scaling method (see @TrajScale::ScaleMethod enum)
// The scaling is returned as matrix T_sc that premultiplies a position to
// scale it, i.e. Y_scaled = T_sc * Y
// The scaling is also invertible, i.e. Y = inv(T_sc) * Y_scaled

namespace as64_
{

namespace gmp_
{

class TrajScale
{

// ===================================
// =======  Public Functions  ========
// ===================================
public:

  /** Constructor.
   * @param[in] n_dof: degrees of freedom.
   */
  TrajScale(unsigned n_dof)
  {
    this->Y0d = arma::vec().zeros(n_dof);
    this->Ygd = arma::vec().ones(n_dof);
    //Initialization required for 'setScaleMethod'
    //will be overwirtten then on the 'setInitFinalPos' call
    this->Y0 = arma::vec().zeros(n_dof);
    this->Yg = arma::vec().ones(n_dof);

    this->setNominalStartFinalPos(this->Y0d, this->Ygd);
    this->setNewStartFinalPos(this->Y0, this->Yg);
  }

  /** Returns the matrix that scales to the new init/target position, according to the set scaling method.
   * @return: The scaling matrix to the new init/target position.
   */
  arma::mat getScaling()
  {
    return this->T_sc;
  }

  /** Returns the inverse of the matrix that scales to the new init/target position, according to the set scaling method.
   * @return: The inverse scaling matrix to the new init/target position.
   */
  arma::mat getInvScaling()
  {
    return this->inv_T_sc;
  }

  /** Set the new start and final position.
   * @param[in] Y0: new start position.
   * @param[in] Yg: new final position.
   */
  void setNewStartFinalPos(const arma::vec &Y0, const arma::vec &Yg)
  {
    this->Y0 = Y0;
    this->Yg = Yg;

    this->T_sc = this->calcScaling();
    this->inv_T_sc = this->calcInvScaling();
  }

  /** Set the nominal start and final position.
   * @param[in] Y0d: nominal start position.
   * @param[in] Ygd: nominal final position.
   */
  void setNominalStartFinalPos(const arma::vec &Y0d, const arma::vec &Ygd)
  {
    this->Y0d = Y0d;
    this->Ygd = Ygd;

    this->T_sc = this->calcScaling();
    this->inv_T_sc = this->calcInvScaling();
  }

  /** Returns the number of DoFs.
   * @return: number of DoFs.
   */
  unsigned numOfDoFs()
  {
    return this->Y0d.size();
  }

protected:

  // =====================================
  // =======  Abstract Functions  ========
  // =====================================

  virtual arma::mat calcScaling() const = 0;
  virtual arma::mat calcInvScaling() const = 0;

  // =======================================
  // =======  Protected Properties  ========
  // =======================================

  arma::vec Y0d; ///< Nominal start position.
  arma::vec Ygd; ///< Nominal final position.
  arma::vec Y0; ///< New start position.
  arma::vec Yg; ///< New final position.

  arma::mat T_sc; //scaling matrix, calculated based on the set scale method
  arma::mat inv_T_sc; //inverse scaling matrix, calculated based on the set scale method

}; // TrajScale

} // namespace gmp_

} // namespace as64_

#endif // GMP_TRAJ_SCALE_H
