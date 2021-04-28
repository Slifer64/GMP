// N-DoF GMP Update class
//

class GMP_nDoF_Update
{

public:

  /**GMP constructor.
   *  @param[in] gmp: n_DoF dmp.
   */
  GMP_nDoF_Update(const std::shared_ptr<gmp_::GMP_nDoF> &gmp)
  {
      this->gmp = gmp;

      this->initSigmaw();
      this->rv = 1.0;

      this->enableSigmawUpdate(false);
  }

  void initSigmaw()
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

  void initSigmaWfromMsr(const arma::rowvec &x_data)
  {
      unsigned n_data = x_data.size();
      arma::mat H(gmp->numOfKernels(), n_data);
      for (int j=0; j<n_data; j++) H.col(j) = gmp->regressVec(x_data(j));

      this->Sigma_w = arma::inv(H*H.t());
  }

  void enableSigmawUpdate(flag)
  {
      this->enable_Sigma_w_update = flag;
  }

  void setSigmaW(const arma::mat &Sw)
  {
    this->Sigma_w = Sw;
  }

  arma::mat getSigmaW() const
  {
    return this->Sigma_w;
  }

  double setMsrNoiseVar(double rv)
  {
    this->rv = rv;
  }

  // ==================================================
  // ===============   Online update  =================
  // ==================================================

  void updatePos(double x, const arma::vec &y, double r_n=-1)
  {
    if (r_n < 0) r_n=this->rv;
    this->updateWeights(gmp_::Phase(x,0,0), y, gmp_::UPDATE_TYPE.POS, r_n);

  }

  void updateVel(double x, double x_dot, const arma::vec &y_dot, double r_n=-1)
  {
    if (r_n < 0) r_n=this->rv;
    this->updateWeights(gmp_::Phase(x,x_dot,0), y_dot, gmp_::UPDATE_TYPE.VEL, r_n);

  }

  void updateAccel(double x, double x_dot, double x_ddot, const arma::vec &y_ddot, r_n=-1)
  {
    if (r_n < 0) r_n=this->rv;
    this->updateWeights(gmp_::Phase(x,x_dot,x_ddot), y_ddot, gmp_::UPDATE_TYPE.ACCEL, r_n);
  }


  /** Updates the weights based on 'n' measurements.
   *  @param[in] s: 1 x n vector of type gmp_::Phase, where the j-th entry is the phase for the j-th measurement.
   *  @param[in] Z: n_dof x n matrix, where the j-th column is the j-th measurement.
   *  @param[in] type: 1 x n vector of type gmp_::UPDATE_TYPE, where the j-th entry has the update type for the j-th measurement.
   *  @param[in] r_n: 1 x n vector, where the j-th entry is the noise variance for the j-th measurement (optional, default=this->rv).
   */
  void updateWeights(const std::vector<gmp_::Phase> &s, const arma::mat &Z, const std::vector<gmp_::UPDATE_TYPE> &type, arma::rowvec r_n={})
  {
    if (r_n.is_empty()) r_n = {this->rv};

    unsigned n = Z.n_cols;
    unsigned n_ker = gmp->numOfKernels();
    arma::mat sc = gmp->getScaling();
    arma::mat inv_sc = gmp->getInvScaling();

    if (r_n.size() == 1) r_n = arma::repmat(r_n(0), 1,n);
    arma::mat Rn = arma::diagmat(r_n);

    arma::mat H(n_ker, n);

    for (int j=0; j<n; j++)
    {
        if (type(j) == gmp_::UPDATE_TYPE.POS)
            H.col(j) = gmp->regressVec(s(j).x);
            Z.col(j) = Z.col(j) - gmp->Y0 + sc*gmp->Y0d;
        else if (type(j) == gmp_::UPDATE_TYPE.VEL)
            H.col(j) = gmp->regressVecDot(s(j).x, s(j).x_dot);
        else % (type == gmp_::UPDATE_TYPE.ACCEL)
            H.col(j) = gmp->regressVecDDot(s(j).x, s(j).x_dot, s(j).x_ddot);
        }
        //Z.col(j) = Z.col(j)/sc;
    }
    Z = inv_sc*Z;

    C = this->Sigma_w*H;
    B = arma::inv_sympd(H.t()*this->Sigma_w*H + Rn) * C.t();
    gmp->W = gmp->W + (Z - gmp->W*H)*B;

    if (this->enable_Sigma_w_update) this->Sigma_w = this->Sigma_w - C*B;
  }

private:

  std::shared_ptr<gmp_::GMP_nDoF> gmp; ///< n_DoF GMP

  arma::mat Sigma_w; ///< weights covariance for trajectory update
  double rv; ///< noise variance for trajectory update
  bool enable_Sigma_w_update;

}
