#include <math_lib/svf.h>

namespace as64_
{

namespace math_
{

/**
* @brief SingularValueFilter::init
* @param sigma_min, minimum allowable eigen value.
* @param shape_f, the greater the shape factor the closer are the filtered eigenvalues to the initial ones.
*/
SingularValueFilter::SingularValueFilter(double sigma_min, double shape_f)
{
  set_sigma_min(sigma_min);
  set_shape_factor(shape_f);
}


/**
* @brief SingularValueFilter::set_sigma_min
* @param sigma_min, minimum allowable eigen value.
*/
void SingularValueFilter::set_sigma_min(double sigma_min)
{
  sigma0 = sigma_min;
}


/**
* @brief SingularValueFilter::set_shape_factor
* @param shape_f, the greater the shape factor the closer are the filtered eigenvalues to the initial ones.
*/
void SingularValueFilter::set_shape_factor(double shape_f)
{
  v = shape_f;
}

double SingularValueFilter::get_sigma_min() const
{
  return sigma0;
}

double SingularValueFilter::get_shape_factor() const
{
  return v;
}

double SingularValueFilter::filter_eig_val(double sigma) const
{
  double filt_sigma = ( std::pow(sigma,3.0) + v*std::pow(sigma,2.0) + 2*sigma + 2*sigma0 ) / (std::pow(sigma,2.0) + v*sigma + 2);

  return filt_sigma;
}

Eigen::MatrixXd SingularValueFilter::inv(Eigen::MatrixXd M) const
{
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(M, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::MatrixXd S = svd.singularValues().asDiagonal();
  Eigen::MatrixXd U = svd.matrixU();
  Eigen::MatrixXd V = svd.matrixV();

  // Eigen::MatrixXd M_reconstr = U*S*V.transpose();

  for (int i=0;i<S.cols(); i++) S(i,i) = 1/filter_eig_val(S(i,i));

  return V*S*U.transpose();
}

} // namespace math_

} // namespace as64_
