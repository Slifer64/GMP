#ifndef GMP_TRAJ_ROT_MIN_SCALE_H
#define GMP_TRAJ_ROT_MIN_SCALE_H

#include <gmp_lib/TrajScale/TrajScale.h>
#include <gmp_lib/math/math.h>

namespace as64_
{

namespace gmp_
{

class TrajScale_Rot_min: public TrajScale
{
public:

  /** Constructor.
   *  @param[in] n_dof: degrees of freedom.
   */
  TrajScale_Rot_min(): TrajScale(3) {}

protected:

  arma::mat calcScaling() const
  {
    arma::vec nd = this->Ygd - this->Y0d;  nd = nd/arma::norm(nd);
    arma::vec n = this->Yg - this->Y0;  n = n/arma::norm(n);
    double dot_n_nd = arma::dot(n,nd);
    arma::mat R;
    if (std::fabs(std::fabs(dot_n_nd) - 1) < 1e-14)
      R = arma::mat().eye(3,3);
    else
    {
      arma::vec k = arma::cross(nd,n);
      double theta = std::acos(dot_n_nd);
      R = gmp_::axang2rotm( arma::join_vert(k, arma::vec({theta})) );
    }

    return R*( arma::norm(this->Yg - this->Y0)/arma::norm(this->Ygd - this->Y0d) );
  }

  arma::mat calcInvScaling() const
  {
    arma::vec nd = this->Ygd - this->Y0d;  nd = nd/arma::norm(nd);
    arma::vec n = this->Yg - this->Y0;  n = n/arma::norm(n);
    double dot_n_nd = arma::dot(n,nd);
    arma::mat R;
    if (std::fabs(std::fabs(dot_n_nd) - 1) < 1e-14)
      R = arma::mat().eye(3,3);
    else
    {
      arma::vec k = arma::cross(nd,n);
      double theta = std::acos(dot_n_nd);
      R = gmp_::axang2rotm( arma::join_vert(k, arma::vec({theta})) );
    }
    return R.t()*( arma::norm(this->Ygd - this->Y0d)/arma::norm(this->Yg - this->Y0) );
  }

}; // TrajScale_Rot_min

} // namespace gmp_

} // namespace as64_

#endif // GMP_TRAJ_ROT_MIN_SCALE_H
