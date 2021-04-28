#ifndef DYNAMIC_TIME_WARPING_64_H
#define DYNAMIC_TIME_WARPING_64_H

#include <cstring>
#include <cmath>
#include <exception>
#include <algorithm>

#include <armadillo>

namespace as64_
{

namespace math_
{

typedef double (*distance_function)(const arma::vec &x1, const arma::vec &x2);

/** function Dynamic Time Warping
 *  Performs dynamic time warping on two signals with dimensions D x N1 and
 *  D x N2 respectively. D is the dimensionality of the input data and N1 and
 *  N2 the number of points in each signal.
 *  @param[in] s: D x N1 matrix with the first input signal.
 *  @param[in] t: D x N2 matrix with the second input signal.
 *  @param[in] w: Window to search for matching points in DTW (optional, default = Inf).
 *  @param[in] dist_fun_h: Pointer to the distance function for the input data (optional, default = 'norm').
 *  @param[out] d: The distance found by the DTW.
 *  @param[out] ind_s: Indices of the first signal after the DTW.
 *  @param[out] ind_t: Indices of the first signal after the DTW.
 */
double dtw(const arma::mat &s, const arma::mat &t, arma::rowvec &ind_s, arma::rowvec &ind_t, int w = -1, distance_function dist_fun_h=arma::norm);

} // namespace math_

} // namespace as64_

#endif // DYNAMIC_TIME_WARPING_64_H
