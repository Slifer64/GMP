#ifndef AS64_gmp_MATH_LIB_H
#define AS64_gmp_MATH_LIB_H

#include <Eigen/Dense>
#include <armadillo>

namespace as64_
{

namespace gmp_
{

Eigen::Vector4d quatInv(const Eigen::Vector4d &quat);
arma::vec quatInv(const arma::vec &quat);

Eigen::Matrix4d quat2mat(const Eigen::Vector4d &quat);

Eigen::Vector4d quatProd(const Eigen::Vector4d &quat1, const Eigen::Vector4d &quat2);
arma::vec quatProd(const arma::vec &quat1, const arma::vec &quat2);

Eigen::Vector4d quatExp(const Eigen::Vector3d &v_rot, double zero_tol=1e-16);
arma::vec quatExp(const arma::vec &v_rot, double zero_tol=1e-16);

Eigen::Vector3d quatLog(const Eigen::Vector4d &quat, double zero_tol=1e-16);
arma::vec quatLog(const arma::vec &quat, double zero_tol=1e-16);

} // namespace gmp_

} // namespace as64_

#endif // AS64_gmp_MATH_LIB_H
