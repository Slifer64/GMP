#ifndef SINGULAR_VALUE_FILTER_H
#define SINGULAR_VALUE_FILTER_H

#include <cmath>
#include <Eigen/Dense>

namespace as64_
{

namespace math_
{

class SingularValueFilter
{
public:
    SingularValueFilter(double sigma_min = 1e-3, double shape_f = 5);

    void set_sigma_min(double sigma_min);
    void set_shape_factor(double shape_f);

    double get_sigma_min() const;
    double get_shape_factor() const;

    double filter_eig_val(double sigma) const;

    Eigen::MatrixXd inv(Eigen::MatrixXd M) const;
private:
    double sigma0; // minimum allowable eigen value
    double v; // shape factor
};

} // namespace math_

} // namespace as64_

#endif // SINGULAR_VALUE_FILTER_H
