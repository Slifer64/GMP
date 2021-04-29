#ifndef GMP_TRAJ_ROT_WB_SCALE_H
#define GMP_TRAJ_ROT_WB_SCALE_H

#include <gmp_lib/TrajScale/TrajScale.h>
#include <gmp_lib/math/math.h>

namespace as64_
{

namespace gmp_
{

class TrajScale_Rot_wb : public TrajScale
{
public:
  /** Constructor.
   *  @param[in] n_dof: degrees of freedom.
   */
  TrajScale_Rot_wb(): TrajScale(3)
  {
    this->setNominalStartFinalPos(this->Y0d, this->Ygd);
    this->setNewStartFinalPos(this->Y0, this->Yg);
    this->setWorkBenchNormal({0,0,1});
  }

  /** Sets the normal to the work bench vector. Applies only to 'ROT_WB_SCALE'.
   *  @param[in] n_wb: 3x1 vector of the normal to the workbench.
   */
  void setWorkBenchNormal(const arma::vec &n_wb)
  {
    this->n_wb = n_wb;

    this->calcScaling();
    this->calcInvScaling();
  }

  ScaleType getScaleType() const { return TrajScale::ScaleType::ROT_WB_SCALE; }

protected:

  arma::mat calcScaling() const
  {
    arma::vec nd = this->Ygd - this->Y0d;  nd = nd/arma::norm(nd);
    arma::vec n = this->Yg - this->Y0;  n = n/arma::norm(n);
    double dot_n_nd = arma::dot(n,nd);
    arma::mat R;
    if (std::fabs(std::fabs(dot_n_nd) - 1) < 1e-12)
      R = arma::mat().eye(3,3);
    else
    {
      arma::vec p = (n - nd) / arma::norm(n - nd);
      arma::vec k = this->n_wb - arma::dot(p,this->n_wb)*p;
      k = k / arma::norm(k);

      arma::vec nd2 = nd - arma::dot(k,nd)*k;
      arma::vec n2 = n - arma::dot(k,n)*k;
      double theta = std::acos( arma::dot(n2,nd2) / (arma::norm(n2)*arma::norm(nd2)) );

      R = gmp_::axang2rotm( arma::join_vert(k, arma::vec({theta})) );
    }

    return R * ( arma::norm(this->Yg - this->Y0)/arma::norm(this->Ygd - this->Y0d) );

  }

  arma::mat calcInvScaling() const
  {
    return arma::inv(this->calcScaling());
  }


private:
  arma::vec n_wb; ///< workbench normal

}; // TrajScale_Rot_wb

} // namespace gmp_

} // namespace as64_

#endif // GMP_TRAJ_ROT_WB_SCALE_H
