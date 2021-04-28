#include <gmp_lib/GMP/GMPo.h>
#include <gmp_lib/math/quaternions.h>

namespace as64_
{

namespace gmp_
{

double GMPo::zero_tol = 1e-7;

// ===================================
// =======  Public Functions  ========
// ===================================

GMPo::GMPo(arma::uvec N_kernels, arma::vec D, arma::vec K, double kernels_std_scaling):
          GMP_nDoF(3, N_kernels, D, K, kernels_std_scaling)
{
  #ifdef GMPo_DEBUG_
  try{
  #endif

  this->setQ0(arma::vec({1, 0, 0, 0}));

  #ifdef GMPo_DEBUG_
  }catch(std::exception &e) { throw std::runtime_error(std::string("[GMPo::GMPo]: ") + e.what()); }
  #endif
}


void GMPo::train(const std::string &train_method, const arma::rowvec &Time, const arma::mat &Quat_data, arma::vec *train_error)
{
  #ifdef GMPo_DEBUG_
  try{
  #endif

  unsigned n_data = Time.size();
  this->setQ0(Quat_data.col(0));
  Qd0 = Quat_data.col(0);

  arma::mat qd_data(3, n_data);
  for (int j=0; j<n_data; j++)
  {
    arma::vec Q1 = GMPo::quatTf(Quat_data.col(j), this->Q0);
    qd_data.col(j) = quatLog(Q1);
  }

  if (train_error) GMP_nDoF::train(train_method, Time, qd_data, train_error);
  else GMP_nDoF::train(train_method, Time, qd_data);

  #ifdef GMPo_DEBUG_
  }catch(std::exception &e) { throw std::runtime_error(std::string("[GMPo::train]: ") + e.what()); }
  #endif
}


void GMPo::setQ0(const arma::vec &Q0)
{
  #ifdef GMPo_DEBUG_
  try{
  #endif

  this->Q0 = Q0;

  #ifdef GMPo_DEBUG_
  }catch(std::exception &e) { throw std::runtime_error(std::string("[GMPo::setY0]: ") + e.what()); }
  #endif
}


void GMPo::setQg(const arma::vec &Qg)
{
  #ifdef GMPo_DEBUG_
  try{
  #endif

  arma::vec qg = GMPo::quat2q(Qg, this->Q0);
  GMP_nDoF::setGoal(qg);

  #ifdef GMPo_DEBUG_
  }catch(std::exception &e) { throw std::runtime_error(std::string("[GMPo::setGoal]: ") + e.what()); }
  #endif
}


arma::vec GMPo::getRotVel(const arma::vec &Q) const
{
  arma::vec Q1 = GMPo::quatTf(Q, this->Q0);
  arma::vec qdot = this->getYdot();
  return GMPo::qdot2rotVel(qdot, Q1);
}

arma::vec GMPo::getRotAccel(const arma::vec &Q, double tau_dot, const arma::vec &Yc_dot) const
{
  arma::vec Q1 = GMPo::quatTf(Q, this->Q0);
  arma::vec qddot = this->getYddot(Yc_dot);
  arma::vec rotVel = this->getRotVel(Q);
  return GMPo::qddot2rotAccel(qddot, rotVel, Q1);
}

arma::vec GMPo::calcRotAccel(const gmp_::Phase &s, const arma::vec &Q, const arma::vec &rotVel, const arma::vec &Qg,
                               const arma::vec &Yc, const arma::vec &Zc, const arma::vec &Yc_dot) const
{
  arma::vec D = arma::vec({this->gmp[0]->D, this->gmp[1]->D, this->gmp[2]->D});
  arma::vec K = arma::vec({this->gmp[0]->K, this->gmp[1]->K, this->gmp[2]->K});

  arma::vec qg = GMPo::quat2q(Qg, this->Q0);

  arma::vec Q1 = GMPo::quatTf(Q, this->Q0);
  arma::vec q = GMPo::quat2q(Q, this->Q0);
  arma::vec invQ1 = quatInv(Q1);

  arma::mat JQq = GMPo::jacobQq(Q1);
  arma::mat JQq_dot = GMPo::jacobDotQq(Q1, rotVel);

  arma::vec fo(3);
  for (int i=0;i<3;i++) fo(i) = this->gmp[i]->shapeAttractor(s);

  arma::vec qdot = GMPo::rotVel2qdot(rotVel, Q1);
  arma::vec qddot = K%(qg - q) - D%(qdot + Yc) + fo + Yc_dot + Zc;

  arma::vec rotAccel = 2*quatProd(JQq_dot*qdot + JQq*qddot, invQ1);

  return rotAccel.subvec(1,3);
}

arma::vec GMPo::getY(const arma::vec &Q) const
{
  return GMPo::quat2q(Q, this->Q0);
}

arma::vec GMPo::getZ(const arma::vec &rotVel, const arma::vec &Q) const
{
  arma::vec Q1 = GMPo::quatTf(Q, this->Q0);
  arma::vec qdot = GMPo::rotVel2qdot(rotVel, Q1);
  arma::vec z = qdot;

  return z;
}


void GMPo::exportToFile(const std::string &filename) const
{
  FileIO fid(filename, FileIO::out | FileIO::trunc);
  writeToFile(fid, "");
}

std::shared_ptr<GMPo> GMPo::importFromFile(const std::string &filename)
{
  std::shared_ptr<GMPo> gmp( new gmp_::GMPo(arma::uvec({2}), arma::vec({1}), arma::vec({1}), 1) );
  FileIO fid(filename, FileIO::in);
  gmp->readFromFile(fid, "");
  return gmp;
}

void GMPo::writeToFile(FileIO &fid, const std::string &prefix) const
{
  GMP_nDoF::writeToFile(fid, prefix);
  fid.write(prefix+"Qd0", this->Qd0);
}

void GMPo::readFromFile(FileIO &fid, const std::string &prefix)
{
  GMP_nDoF::readFromFile(fid, prefix);
  fid.read(prefix+"Qd0", this->Qd0);
  this->Q0 = this->Qd0;
}


arma::vec GMPo::getQd(double x) const
{
  return gmp_::GMPo::q2quat(getYd(x), this->Q0);
}

arma::vec GMPo::getRotVeld(double x, double x_dot) const
{
  arma::vec q_dot = getYdDot(x, x_dot);
  arma::vec Q1 = gmp_::GMPo::quatTf(getQd(x),Qd0);
  return gmp_::GMPo::qdot2rotVel( q_dot, Q1 );
}

arma::vec GMPo::getRotAcceld(double x, double x_dot, double x_ddot) const
{
  arma::vec q_ddot = getYdDDot(x, x_dot, x_ddot);
  arma::vec Q1 = gmp_::GMPo::quatTf(getQd(x),Qd0);
  arma::vec rotVel = getRotVeld(x, x_dot);
  return gmp_::GMPo::rotAccel2qddot(q_ddot, rotVel, Q1);
}


// ====================================================
// ****************************************************
// ************      Static Functions      ************
// ****************************************************
// ====================================================

arma::vec GMPo::quatTf(const arma::vec &Q, const arma::vec &Q0)
{
  return quatProd(Q, quatInv(Q0));
}


arma::vec GMPo::quat2q(const arma::vec &Q, const arma::vec &Q0)
{
  return quatLog(GMPo::quatTf(Q,Q0));
}


arma::vec GMPo::q2quat(const arma::vec &q, const arma::vec &Q0)
{
  return quatProd( quatExp(q), Q0);
}


arma::vec GMPo::rotVel2qdot(const arma::vec &rotVel, const arma::vec &Q1)
{
  arma::mat JqQ = GMPo::jacobqQ(Q1);
  return 0.5 * JqQ * quatProd( arma::join_vert(arma::vec({0}), rotVel), Q1 );
}


arma::vec GMPo::qdot2rotVel(const arma::vec &qdot, const arma::vec &Q1)
{
  arma::mat JQq = GMPo::jacobQq(Q1);
  arma::vec rotVel = 2 * quatProd( JQq*qdot, quatInv(Q1) );
  return rotVel.subvec(1,3);
}


arma::vec GMPo::rotAccel2qddot(const arma::vec &rotAccel, const arma::vec &rotVel, const arma::vec &Q1)
{
  arma::vec rotVelQ = arma::join_vert(arma::vec({0}), rotVel);
  arma::vec rotAccelQ = arma::join_vert(arma::vec({0}), rotAccel);

  arma::mat J = GMPo::jacobqQ(Q1);
  arma::mat Jdot = GMPo::jacobDotqQ(Q1, rotVel);

  return 0.5*(Jdot * quatProd(rotVelQ, Q1) + J * quatProd( rotAccelQ+0.5*quatProd(rotVelQ,rotVelQ), Q1 ) );
}


arma::vec GMPo::qddot2rotAccel(const arma::vec &qddot, const arma::vec &rotVel, const arma::vec &Q1)
{
  arma::vec qdot = GMPo::rotVel2qdot(rotVel, Q1);
  arma::vec invQ1 = quatInv(Q1);
  arma::mat J = GMPo::jacobQq(Q1);
  arma::mat Jdot = GMPo::jacobDotQq(Q1, rotVel);

  arma::vec rotAccel = 2 * ( quatProd( Jdot*qdot+J*qddot, invQ1 ) );
  return rotAccel.subvec(1,3);
}


arma::mat GMPo::jacobQq(const arma::vec &Q1)
{
  arma::mat JQq(4,3);

  if ( (1-std::fabs(Q1(0))) <= GMPo::zero_tol)
  {
    JQq.row(0) = arma::rowvec().zeros(3);
    JQq.submat(1,0,3,2) = arma::mat().eye(3,3);
    return JQq;
  }

  double w = Q1(0);
  arma::vec v = Q1.subvec(1,3);
  double norm_v = arma::norm(v);
  arma::vec eta = v / norm_v;
  double s_th = norm_v;
  double c_th = w;
  double th = std::atan2(s_th, c_th);
  arma::mat Eta = eta*eta.t();

  JQq.row(0) = -0.5 * s_th * eta.t();
  JQq.submat(1,0,3,2) = 0.5 * ( (arma::mat().eye(3,3) - Eta)*s_th/th + c_th*Eta );
  return JQq;
}


arma::mat GMPo::jacobqQ(const arma::vec &Q1)
{
  arma::mat JqQ(3,4);

  if ( (1-std::fabs(Q1(0))) <= GMPo::zero_tol)
  {
    JqQ.col(0) = arma::vec().zeros(3);
    JqQ.submat(0,1,2,3) = arma::mat().eye(3,3);
    return JqQ;
  }

  double w = Q1(0);
  arma::vec v = Q1.subvec(1,3);
  double norm_v = arma::norm(v);
  arma::vec eta = v / norm_v;
  double s_th = norm_v;
  double c_th = w;
  double th = std::atan2(s_th, c_th);

  JqQ.col(0) = 2*eta*(th*c_th - s_th)/std::pow(s_th,2);
  JqQ.submat(0,1,2,3) = 2*arma::mat().eye(3,3)*th/s_th;
  return JqQ;
}


arma::mat GMPo::jacobDotqQ(const arma::vec &Q1, const arma::vec &rotVel)
{
  arma::mat JqQ_dot(3,4);

  arma::vec qdot = GMPo::rotVel2qdot(rotVel, Q1);

  if ( (1-std::fabs(Q1(0))) <= GMPo::zero_tol)
  {
    JqQ_dot.col(0) = -qdot/3;
    JqQ_dot.submat(0,1,2,3) = arma::mat().zeros(3,3);
    return JqQ_dot;
  }

  double w = Q1(0);
  arma::vec v = Q1.subvec(1,3);
  double norm_v = arma::norm(v);
  arma::vec eta = v / norm_v;
  double s_th = norm_v;
  double c_th = w;
  double th = std::atan2(s_th, c_th);
  arma::mat Eta = eta*eta.t();
  double temp = (th*c_th-s_th)/std::pow(s_th,2);

  JqQ_dot.col(0) = ((-th/s_th - 2*c_th*temp/s_th)*Eta + temp*(arma::mat().eye(3,3)-Eta)/th)*qdot;
  JqQ_dot.submat(0,1,2,3) = (-temp*arma::dot(eta,qdot))*arma::mat().eye(3,3);
  return JqQ_dot;
}


arma::mat GMPo::jacobDotQq(const arma::vec &Q1, const arma::vec &rotVel)
{
  arma::mat JQq_dot(4,3);

  arma::vec qdot = GMPo::rotVel2qdot(rotVel, Q1);

  if ( (1-std::fabs(Q1(0))) <= GMPo::zero_tol)
  {
    JQq_dot.row(0) = -qdot.t()/4;
    JQq_dot.submat(1,0,3,2) = arma::mat().zeros(3,3);
    return JQq_dot;
  }

  double w = Q1(0);
  arma::vec v = Q1.subvec(1,3);
  double norm_v = arma::norm(v);
  arma::vec eta = v / norm_v;
  double s_th = norm_v;
  double c_th = w;
  double th = std::atan2(s_th, c_th);
  arma::mat Eta = eta*eta.t();
  arma::mat I_eta = arma::mat().eye(3,3) - Eta;
  double temp = ((th*c_th-s_th)/std::pow(th,2));

  JQq_dot.row(0) = -0.25 * qdot.t() * (c_th*Eta + (s_th/th)*I_eta);
  JQq_dot.submat(1,0,3,2) = (0.25*arma::dot(eta,qdot))*( temp*I_eta - s_th*Eta ) + 0.25*temp*( eta*(qdot.t()*I_eta) + (I_eta*qdot)*eta.t() );

  return JQq_dot;
}


} // namespace gmp_

} // namespace as64_
