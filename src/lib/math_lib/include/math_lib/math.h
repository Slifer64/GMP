#ifndef MATH_LIB_MATH__64_H
#define MATH_LIB_MATH__64_H

#include <Eigen/Dense>
#include <armadillo>

namespace as64_
{

namespace math_
{

Eigen::Matrix3d vec2ssMat(const Eigen::Vector3d &v);
arma::mat vec2ssMat(const arma::vec &v);

Eigen::Vector4d rotm2quat(Eigen::Matrix3d rotm, bool is_rotm_orthonormal=true);
arma::vec rotm2quat(const arma::mat &rotm, bool is_rotm_orthonormal=true);
Eigen::Vector4d mat2quat(Eigen::Matrix3d rotm);

Eigen::Matrix3d quat2rotm(Eigen::Vector4d quat);
arma::mat quat2rotm(const arma::vec &quat);

Eigen::Vector4d rotm2axang(Eigen::Matrix3d rotm);
arma::vec rotm2axang(const arma::mat &rotm);

Eigen::Matrix3d axang2rotm(Eigen::Vector4d axang);
arma::mat axang2rotm(const arma::vec &axang);

Eigen::Vector4d axang2quat(Eigen::Vector4d axang);
arma::vec axang2quat(const arma::vec &axang);

Eigen::Vector4d quat2axang(Eigen::Vector4d quat);
arma::vec quat2axang(const arma::vec &quat);

Eigen::MatrixXd inv(const Eigen::MatrixXd &M);
arma::mat inv(const arma::mat &M);

/**
 * \brief Returns a quaternion that represents the same rotation as "quat1" and as
 *  a 4D vector creates the smaller angle with "quat2"
 */
Eigen::Vector4d getClosestQuat(const Eigen::Vector4d &quat1, const Eigen::Vector4d &quat2);
arma::vec getClosestQuat(const arma::vec &quat1, const arma::vec &quat2);

/** \brief Returns the zero, first and second order derivatives of a fifth order polynomial.
 * @param[in] t Timestamp.
 * @param[in] p0 Initial state.
 * @param[in] pT Final state.
 * @param[in] totalTime Total time duration.
 * @return Matrix with the zero, first and second order derivatives in each row for timestamp 't'.
 */
arma::mat get5thOrder(double t, arma::vec p0, arma::vec pT, double totalTime);

} // namespace math_

} // namespace as64_

#endif // MATH_LIB_MATH__64_H
