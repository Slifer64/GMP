#include <gmp_lib/GMP/GMP_Update.h>

namespace as64_
{

namespace gmp_
{

  GMP_Update::GMP_Update(gmp_::GMP *gmp)
  {
      this->gmp = gmp;

      this->initSigmaw();
      this->rv = 1.0;

      this->enableSigmawUpdate(false);
  }

  void GMP_Update::initSigmaw()
  {
      unsigned N_kernels = gmp->numOfKernels();
      // arma::mat S = arma::mat().zeros(N_kernels,N_kernels);
      // for (int i=0; i<N_kernels; i++)
      // {
      //     for (int j=i+1; j<N_kernels; j++)
      //     {
      //         S(i,j) = std::exp(-0.2 * std::abs(i-j));
      //     }
      // }
      // S = S + S.t() + arma::mat().eye(N_kernels,N_kernels);
      // this->Sigma_w = S;
      this->Sigma_w = arma::mat().eye(N_kernels, N_kernels);
  }

  void GMP_Update::initExpSigmaw(double decay_rate)
  {      
    unsigned N_kernels = this->gmp->numOfKernels();
    arma::mat S = arma::mat().zeros(N_kernels, N_kernels);
    for (int i=0; i<N_kernels; i++)
    {
      for (int j=i+1; j<N_kernels; j++) S(i,j) = std::exp(-decay_rate * std::abs(i-j));
    }
    this->Sigma_w = S + S.t() + arma::mat().eye(N_kernels, N_kernels); 
  }

  void GMP_Update::initSigmaWfromMsr(const arma::rowvec &x_data)
  {
      unsigned n_data = x_data.size();
      arma::mat H(gmp->numOfKernels(), n_data);
      for (int j=0; j<n_data; j++) H.col(j) = gmp->regressVec(x_data(j));

      this->Sigma_w = arma::inv(H*H.t());
  }

  void GMP_Update::enableSigmawUpdate(bool flag)
  {
      this->enable_Sigma_w_update = flag;
  }

  void GMP_Update::setSigmaW(const arma::mat &Sw)
  {
    this->Sigma_w = Sw;
  }

  arma::mat GMP_Update::getSigmaW() const
  {
    return this->Sigma_w;
  }

  void GMP_Update::setMsrNoiseVar(double rv)
  {
    this->rv = rv;
  }

  // ==================================================
  // ===============   Online update  =================
  // ==================================================

  void GMP_Update::updatePos(double x, const arma::vec &y, double r_n)
  {
    if (r_n < 0) r_n=this->rv;
    this->updateWeights({gmp_::Phase(x,0,0)}, y, {gmp_::UPDATE_TYPE::POS}, arma::rowvec({r_n}) );

  }

  void GMP_Update::updateVel(double x, double x_dot, const arma::vec &y_dot, double r_n)
  {
    if (r_n < 0) r_n=this->rv;
    this->updateWeights({gmp_::Phase(x,x_dot,0)}, y_dot, {gmp_::UPDATE_TYPE::VEL}, arma::rowvec({r_n}) );

  }

  void GMP_Update::updateAccel(double x, double x_dot, double x_ddot, const arma::vec &y_ddot, double r_n)
  {
    if (r_n < 0) r_n=this->rv;
    this->updateWeights({gmp_::Phase(x,x_dot,x_ddot)}, y_ddot, {gmp_::UPDATE_TYPE::ACCEL}, arma::rowvec({r_n}) );
  }

  void GMP_Update::updateWeights(const std::vector<gmp_::Phase> &s, arma::mat Z, const std::vector<gmp_::UPDATE_TYPE> &type, arma::rowvec r_n)
  {
    if (r_n.is_empty()) r_n = {this->rv};

    unsigned n = Z.n_cols;
    unsigned n_ker = gmp->numOfKernels();
    arma::mat sc = gmp->getScaling();
    arma::mat inv_sc = gmp->getInvScaling();

    if (r_n.size() == 1) r_n = r_n(0)*arma::rowvec().ones(n);
    arma::mat Rn = arma::diagmat(r_n);

    arma::mat H(n_ker, n);

    for (int j=0; j<n; j++)
    {
        if (type[j] == gmp_::UPDATE_TYPE::POS)
        {
            H.col(j) = gmp->regressVec(s[j].x);
            Z.col(j) = Z.col(j) - gmp->Y0 + sc*gmp->Y0d;
        }
        else if (type[j] == gmp_::UPDATE_TYPE::VEL)
            H.col(j) = gmp->regressVecDot(s[j].x, s[j].x_dot);
        else // (type == gmp_::UPDATE_TYPE::ACCEL)
            H.col(j) = gmp->regressVecDDot(s[j].x, s[j].x_dot, s[j].x_ddot);
        //Z.col(j) = Z.col(j)/sc;
    }
    Z = inv_sc*Z;

    arma::mat C = this->Sigma_w*H;
    arma::mat B = arma::inv_sympd(H.t()*this->Sigma_w*H + Rn) * C.t();
    gmp->W = gmp->W + (Z - gmp->W*H)*B;

    if (this->enable_Sigma_w_update) this->Sigma_w = this->Sigma_w - C*B;
  }

} // namespace gmp_

} // namespace as64_
