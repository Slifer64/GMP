#include <math_lib/quaternions.h>

namespace as64_
{

namespace math_
{

Eigen::Vector4d quatProd(const Eigen::Vector4d &quat1, const Eigen::Vector4d &quat2)
{
  Eigen::Vector4d quat12;

  // quat12 = quat2qmat(quat1) * quat2;
  double n1 = quat1(0);
  Eigen::Vector3d e1 = quat1.segment(1,3);

  double n2 = quat2(0);
  Eigen::Vector3d e2 = quat2.segment(1,3);

  quat12(0) = n1*n2 - e1.dot(e2);
  quat12.segment(1,3) = n1*e2 + n2*e1 + e1.cross(e2);

  return quat12;
}

arma::vec quatProd(const arma::vec &quat1, const arma::vec &quat2)
{
  arma::vec quat12(4);

  double n1 = quat1(0);
  arma::vec e1 = quat1.subvec(1,3);

  double n2 = quat2(0);
  arma::vec e2 = quat2.subvec(1,3);

  quat12(0) = n1*n2 - arma::dot(e1,e2);
  quat12.subvec(1,3) = n1*e2 + n2*e1 + arma::cross(e1,e2);

  return quat12;
}


Eigen::Vector4d quatInv(const Eigen::Vector4d &quat)
{
  Eigen::Vector4d quatI;

  quatI(0) = quat(0);
  quatI.segment(1,3) = - quat.segment(1,3);

  return quatI;
}

arma::vec quatInv(const arma::vec &quat)
{
  arma::vec quatI(4);

  quatI(0) = quat(0);
  quatI.subvec(1,3) = - quat.subvec(1,3);

  return quatI;
}

Eigen::Vector4d quatExp(const Eigen::Vector3d &v_rot, double zero_tol)
{
  Eigen::Vector4d quat;
  double norm_v_rot = v_rot.norm();
  double theta = norm_v_rot;

  if (norm_v_rot > zero_tol)
  {
    quat(0) = std::cos(theta/2);
    quat.segment(1,3) = std::sin(theta/2)*v_rot/norm_v_rot;
  }
  else{
    quat << 1, 0, 0, 0;
  }

  return quat;
}

arma::vec quatExp(const arma::vec &v_rot, double zero_tol)
{
  arma::vec quat(4);
  double norm_v_rot = arma::norm(v_rot);
  double theta = norm_v_rot;

  if (norm_v_rot > zero_tol)
  {
    quat(0) = std::cos(theta/2);
    quat.subvec(1,3) = std::sin(theta/2)*v_rot/norm_v_rot;
  }
  else{
    quat << 1 << 0 << 0 << 0;
  }

  return quat;
}

Eigen::Vector3d quatLog(const Eigen::Vector4d &quat, double zero_tol)
{
  Eigen::Vector3d v = quat.segment(1,3);
  double u = quat(0);

  Eigen::Vector3d omega;
  double v_norm = v.norm();

  if (v_norm > zero_tol) omega = 2*std::atan2(v_norm,u)*v/v_norm;
  else omega = Eigen::Vector3d::Zero();

  return omega;
}

arma::vec quatLog(const arma::vec &quat, double zero_tol)
{
  arma::vec v = quat.subvec(1,3);
  double u = quat(0);

  arma::vec omega(3);
  double v_norm = arma::norm(v);

  if (v_norm > zero_tol) omega = 2*std::atan2(v_norm,u)*v/v_norm;
  else omega = arma::vec().zeros(3);

  return omega;
}

Eigen::Matrix4d quat2qmat(const Eigen::Vector4d &quat)
{
  Eigen::Matrix4d mat;
  Eigen::Vector3d e = quat.segment(1,3);
  double n = quat(0);

  mat << n, -e.transpose(), e, n*Eigen::Matrix3d::Identity() + vec2ssMat(e);

  return mat;
}

// If quat1 and quat2 were positions, this would perform quat1-quat2.
// The result is the amount of rotation needed to go from quat2 to quat1, i.e. quatD*quat2 = quat1
// Equivalently, the result is quaternion corresponding to the angular velocity which takes us from quat2 to quat1 in unit time.
Eigen::Vector4d quatDiff(const Eigen::Vector4d &quat1, const Eigen::Vector4d &quat2)
{
  Eigen::Vector4d quatD;

  quatD = quatProd(quat1, quatInv(quat2));

  return quatD;
}

arma::vec quatDiff(const arma::vec &quat1, const arma::vec &quat2)
{
  arma::vec quatD;

  quatD = quatProd(quat1, quatInv(quat2));

  return quatD;
}


} // namespace math_

} // namespace as64_
