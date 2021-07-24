#ifndef GMP_TRAJ_PROP_SCALE_H
#define GMP_TRAJ_PROP_SCALE_H

#include <gmp_lib/TrajScale/TrajScale.h>

namespace as64_
{

namespace gmp_
{

class TrajScale_Prop : public TrajScale
{

// ===================================
// =======  Public Functions  ========
// ===================================
public:

  typedef std::shared_ptr<TrajScale_Prop> Ptr;

  /** Constructor.
   *  @param[in] n_dof: degrees of freedom.
   */
  TrajScale_Prop(unsigned n_dof): TrajScale(n_dof)
  {
    this->setNominalStartFinalPos(this->Y0d, this->Ygd);
    this->setNewStartFinalPos(this->Y0, this->Yg);
  }

  enum ScaleType getScaleType() const { return TrajScale::ScaleType::PROP_SCALE; }

protected:

  arma::mat calcScaling() const
  {
    return arma::diagmat( (this->Yg - this->Y0) / (this->Ygd - this->Y0d) );
  }

  arma::mat calcInvScaling() const
  {
    return arma::diagmat( (this->Ygd - this->Y0d) / (this->Yg - this->Y0) );
  }

}; // TrajScale_Prop

} // namespace gmp_

} // namespace as64_

#endif // GMP_TRAJ_PROP_SCALE_H
