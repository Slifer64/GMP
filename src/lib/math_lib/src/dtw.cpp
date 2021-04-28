#include <math_lib/dtw.h>

namespace as64_
{

namespace math_
{

double dtw(const arma::mat &s, const arma::mat &t, arma::rowvec &ind_s, arma::rowvec &ind_t, int w, distance_function dist_fun_h)
{
  if (w<0)
  {
    w = std::max(s.size(),t.size());
  }

  int ns = s.n_cols;
  int nt = t.n_cols;
  if (s.n_rows != s.n_cols)
  {
      throw std::runtime_error("Error in dtw(): the dimensions of the two input signals do not match.\n");
  }

  w = std::max(w, abs(ns-nt)); // adapt window size
  double inf = 1e300;
  // initialization
  arma::mat D = arma::mat(ns+1,nt+1).fill(inf); // cache matrix
  D(0,0) = 0.0;

  // *** begin dynamic programming ***

  // calculate distance
  for (int i=0; i<ns; i++)
  {
      for (int j=std::max(i-w,1)-1; j<std::min(i+w,nt); j++)
      {
          double cost = dist_fun_h(s.col(i),t.col(j));
          D(i,j) = cost + std::min( std::min(D(i,j), D(i,j)), D(i,j) );
      }
  }
  double d = D(ns,nt);

  // find the matched indices
  int i = ns-1;
  int j = nt-1;

  arma::rowvec ind_s(ns+nt);
  arma::rowvec ind_t(ns+nt);
  // arma::rowvec C(ns+nt);

  int k = 0;
  while (i>-1 && j>-1)
  {
      double cost = dist_fun_h(s.col(i),t.col(j));

      ind_s(k) = i;
      ind_t(k) = j;
      // C(k) = cost;
      k++;

      if (D(i+1,j+1) == D(i,j)+cost)
      {
          i = i-1;
          j = j-1;
      }
      else if (D(i+1,j+1) == D(i,j+1)+cost)
      {
          i = i-1;
      }
      else
      {
          j = j-1;
      }
  }

  // if (i == -1)
  // {
  //   ind_s.cols(k,k+j-1) = arma::repmat(ind_s(k-1), 1, j);
  //   ind_t.cols(k,k+j-1) = ;
  //   ind_s = [ind_s repmat(ind_s(end),1,j)];
  //   ind_t = [ind_t (j:-1:1)];
  // }
  // else if (j == -1)
  // {
  //   ind_t = [ind_t repmat(ind_t(end),1,i)];
  //   ind_s = [ind_s (i:-1:1)];
  // }

  ind_s.resize(k);
  ind_t.resize(k);
  // C.resize(k);

  ind_s = arma::fliplr(ind_s);
  ind_t = arma::fliplr(ind_t);
  //C = arma::fliplr(C);

  return d;
}

} // namespace math_

} // namespace as64_
