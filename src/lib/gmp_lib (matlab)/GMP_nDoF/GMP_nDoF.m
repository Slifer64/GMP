#ifndef GMP_LIB_GMP_NDOF_H
#define GMP_LIB_GMP_NDOF_H

// N-DoF GMP class
// Generalized movement primitive.
//

#include <gmp_lib/GMP/GMP_regressor.h>
#include <gmp_lib/TrajScale/TrajScale_Prop.h>
#include <gmp_lib/TrajScale/TrajScale_Rot_min.h>
#include <gmp_lib/TrajScale/TrajScale_Rot_wb.h>

namespace as64_
{

namespace gmp_
{

class GMP_nDoF_Update, GMP_nDoF_IO;

class GMP_nDoF: public GMP_regressor
{

public:

  /** GMP constructor.
   *  @param[in] n_dofs: number of degrees of freedom.
   *  @param[in] N_kernels: the number of kernels
   *  @param[in] kern_std_scale: Scaling for std of kernels (optional, default=1).
   */
  GMP_nDoF(unsigned n_dofs, unsigned N_kernels, double kern_std_scale=1.0) : GMP_regressor(N_kernels, kern_std_scale)
  {
    this->K = 100 * arma::vec().ones(n_dofs);
    this->D = 30 * arma::vec().ones(n_dofs);

    this->W = arma::mat().zeros(n_dofs, N_kernels);

    this->Y0d = arma::vec().zeros(n_dofs);
    this->Ygd = arma::vec().ones(n_dofs);
    this->Y0 = this->Y0d;
    this->Yg = this->Ygd;

    this->setScaleMethod(new TrajScale_Prop(n_dofs));

    this->y_dot = arma::vec().zeros(n_dofs);
    this->z_dot = arma::vec().zeros(n_dofs);
  }


  /** Returns the number of DoFs.
   *  return: number of DoFs.
   */
  unsigned numOfDoFs() const
  {
    return this->W.n_rows;
  }

  /** Returns the number of kernels.
   *  @return: number of kernels.
   */
  unsigned numOfKernels() const
  {
    return this->W.n_cols;
  }

  /** Sets the initial position.
   *  @param[in] y0: initial position.
   */
  void setY0(const arma::vec &y0)
  {
    this->Y0 = y0;
    this->traj_sc->setNewStartFinalPos(this->Y0, this->Yg);
  }

  arma::vec getY0() const
  {
    return this->Y0;
  }

  arma::vec getY0d()
  {
    return this->Y0d;
  }


  /** Set goal position.
   *  @param[in] g: goal position.
   */
  void setGoal(const arma::vec &g)
  {
    this->Yg = g;
    this->traj_sc->setNewStartFinalPos(this->Y0, this->Yg);
  }

  arma::vec getGoal() const
  {
    return this->Yg;
  }

  arma::mat getScaling() const
  {
    return this->traj_sc->getScaling();
  }

  arma::mat getInvScaling() const
  {
    return this->traj_sc->getInvScaling();
  }

  /** Set the trajectory spatial scaling method.
   *  @param[in] traj_scale_obj: pointer to an object of type @TrajScale.
   */
  void setScaleMethod(std::shared_ptr<gmp_::TrajScale> traj_scale_obj)
  {
      if (this->numOfDoFs() != traj_scale_obj->numOfDoFs())
        throw std::runtime_error("[GMP_nDoF::setScaleMethod]: Incompatible number of DoFs...");

      this->traj_sc = traj_scale_obj;
      this->traj_sc->setNominalStartFinalPos(this->Y0, this->Yg);
      this->traj_sc->setNewStartFinalPos(this->Y0, this->Yg);
  }

  /** Trains the GMP.
   * @param[in] train_method: the training method to use, as a string ('LWR', 'LS').
   * @param[in] x: Row vector with the canonical timestamps (in [0 1]) of the training data points.
   * @param[in] yd_data: Matrix with the desired potition for each DoF in each row.
   * @param[out] train_error: The training error expressed as the mse error.
   * \note The timestamps-data need not be sequencial temporarily.
   */
  void train(const std::string &train_method, const arma::rowvec &x, const arma::mat &yd_data, arma::vec *train_err, arma::mat *Sw)
  {
    for (int i=0; i<x.size(); i++)
    { if (x(i)>1 || x(i)<0) throw std::runtime_error("[GMP_nDoF::train]: The training timestamps are not normalized..."); }

    unsigned n_data = x.size();
    unsigned num_ker = this->numOfKernels();
    unsigned n_dofs = this->numOfDoFs();
    unsigned i_end = n_data-1;

    arma::mat H(num_ker, n_data);
    for (j=0; j<n_data; j++) H.col(j) = this->regressVec(x(j));

    if (train_method.compare("LWR") == 0)
        this->W = yd_data*H.t() / arma::repmat(arma::sum(H,1).t(),n_dofs,1);
    else if (train_method.compare("LS") == 0)
        this->W = arma::solve(H.t() yd_data.t()).t(); // yd_data / H;
    else
        throw std::runtime_error("[WSoG::train]: Unsupported training method...");

    this->Y0d = this->W*H.col(0);
    this->Ygd = this->W*H.col(i_end);

    this->setY0(this->Y0d);
    this->setGoal(this->Ygd);

    this->traj_sc->setNominalStartFinalPos(this->Y0d, this->Ygd);
    this->traj_sc->setNewStartFinalPos(this->getY0() this->getGoal());

    if (train_err)
        *train_err = arma::vec().zeros(n_dofs);
        arma::mat err_data = this->W*H - yd_data;
        for (i=0; i<n_dofs; i++) train_err->at(i) = arma::norm(err_data.row(i));
    }

    if (Sw) *Sw = arma::inv(H*H.t());

  }


  /** Calculates the time derivatives of the GMP's states.
   * @param[in] s: Object of type @GMP_phase.
   * @param[in] y: 'y' state of the GMP.
   * @param[in] z: 'z' state of the GMP.
   * @param[in] y_c: coupling term for the dynamical equation of the 'y' state (optional, default=0).
   * @param[in] z_c: coupling term for the dynamical equation of the 'z' state (optional, default=0).
   */
  void update(const gmp_::Phase &s, const arma::vec &y, const arma::vec &z, arma::vec y_c={0}, arma::vec z_c={0})
  {
    unsigned n_dofs = this->numOfDoFs();

    if (y_c.size() == 1) y_c = arma::vec().ones(n_dofs)*y_c(0);
    if (z_c.size() == 1) z_c = arma::vec().ones(n_dofs)*z_c(0);

    arma::vec yd = this->getYd(s.x);
    arma::vec yd_dot = this->getYdDot(s.x, s.x_dot);
    arma::vec yd_ddot = this->getYdDDot(s.x, s.x_dot, s.x_ddot);

    this->y_dot = z + y_c;
    this->z_dot = this->K%(yd - y) + this->D%(yd_dot - z) + yd_ddot + z_c;
  }


  /** Returns the 'y' state time derivative.
   * Call after @update.
   * @return: time derivative of 'y' state.
   */
  arma::vec getYdot() const
  {
    return this->y_dot;
  }


  /** Returns the 'z' state time derivative.
   * Call after @update.
   * @return: time derivative of 'z' state.
   */
  arma::vec getZdot() const
  {
    return this->z_dot;
  }


  /** Returns the GMP's acceleration.
   * Call after @update.
   * @param[in] yc_dot: time derivative of 'y' state coupling (optional, default=0).
   * @return: acceleration.
   */
  arma::vec getYddot(arma::vec yc_dot={0}) const
  {
    unsigned n_dofs = this->numOfDoFs();
    if (yc_dot.size()==1) yc_dot = arma::vec().ones(n_dofs)*yc_dot(0);

    return this->getZdot() + yc_dot;
  }


  /** Calculates the GMP's acceleration.
   * @param[in] s: Object of type @GMP_phase.
   * @param[in] y: 'y' state of the GMP.
   * @param[in] y_dot: time derivative of 'y' state.
   * @param[in] y_c: coupling term for the dynamical equation of the 'y' state (optional, default=0).
   * @param[in] z_c: coupling term for the dynamical equation of the 'z' state (optional, default=0).
   * @param[in] yc_dot: time derivative of the 'y' state coupling (optional, default=0).
   * @return: acceleration.
   */
  arma::vec calcYddot(const gmp_::Phase &s, const arma::vec &y, const arma::vec &y_dot, arma::vec y_c={0}, arma::vec z_c={0}, arma::vec yc_dot={0})
  {
      unsigned n_dofs = this->numOfDoFs();

      if (y_c.size()==1) y_c = arma::vec().ones(n_dofs)*y_c(0);
      if (z_c.size()==1) z_c = arma::vec().ones(n_dofs)*z_c(0);
      if (yc_dot.size()==1) yc_dot = arma::vec().ones(n_dofs)*yc_dot(0);

      yd = this->getYd(s.x);
      yd_dot = this->getYdDot(s.x, s.x_dot);
      yd_ddot = this->getYdDDot(s.x, s.x_dot, s.x_ddot);

      z = y_dot - y_c;
      z_dot = this->K%(yd - y) + this->D%(yd_dot - z) + yd_ddot + z_c;
      return (z_dot + yc_dot);
  }

  /** Returns the scaled desired position.
   * @param[in] x: phase variable.
   * @return: scaled desired position.
   */
  arma::vec getYd(double x) const
  {
    return this->getScaling()*(this->W*this->regressVec(x) - this->Y0d) + this->Y0;
  }


  /** Returns the scaled desired velocity.
   * @param[in] x: phase variable.
   * @param[in] x_dot: 1st time derivative of the phase variable.
   * @return: scaled desired velocity.
   */
  arma::vec getYdDot(double x, double x_dot) const
  {
    return this->getScaling()*this->W*this->regressVecDot(x,x_dot);
  }


  /** Returns the scaled desired acceleration.
   * @param[in] x: phase variable.
   * @param[in] x_dot: 1st time derivative of the phase variable.
   * @param[in] x_ddot: 2nd time derivative of the phase variable.
   * @return: scaled desired acceleration.
   */
  arma::vec getYdDDot(double x, double x_dot, double x_ddot=0) const
  {
    return this->getScaling()*this->W*this->regressVecDDot(x,x_dot,x_ddot);
  }


public: // properties

  std::shared_ptr<gmp_::TrajScale> traj_sc; ///< object of type @TrajScale

  // weights
  arma::mat W; ///< num_DoF x num_Kernels matrix where each row contrains the weights for each DoF

  // impedance params
  arma::vec K; ///< num_DoF x 1 stiffness vector
  arma::vec D; ///< num_DoF x 1 stiffness vector

protected: // properties

friend GMP_nDoF_Update, GMP_nDoF_IO;

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
