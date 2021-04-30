#include <project_name_/utils/utils.h>

arma::wall_clock Timer::timer;


void PRINT_INFO_MSG(const std::string &msg, std::ostream &out)
{
  out << "\033[1m" << "\033[34m" << "[INFO]: " << msg << "\033[0m";
}

void PRINT_CONFIRM_MSG(const std::string &msg, std::ostream &out)
{
  out << "\033[1m" << "\033[32m" << "[INFO]: " << msg << "\033[0m";
}

void PRINT_WARNING_MSG(const std::string &msg, std::ostream &out)
{
  out << "\033[1m" << "\033[33m" << "[WARNING]: " << msg << "\033[0m";
}

void PRINT_ERROR_MSG(const std::string &msg, std::ostream &out)
{
  out << "\033[1m" << "\033[31m" << "[ERROR]: " << msg << "\033[0m";
}

void simulateGMP_nDoF(gmp_::GMP_nDoF::Ptr &gmp, const arma::vec &y0,
                 const arma::vec &yg, double T, double dt,
                 arma::mat &Time, arma::mat &Y_data, arma::mat &dY_data, arma::mat &ddY_data)
{

  int Dim = gmp->numOfDoFs();
  arma::vec y = y0;
  arma::vec dy = arma::vec().zeros(Dim);
  arma::vec ddy = arma::vec().zeros(Dim);
  arma::vec z = arma::vec().zeros(Dim);
  arma::vec dz = arma::vec().zeros(Dim);

  double t = 0.0;
  double t_end = T;
  double tau = t_end;
  double x = 0;
  double x_dot = 1/tau;
  double x_ddot = 0;
  gmp_::Phase s(x, x_dot, x_ddot);

  gmp->setY0(y0);
  gmp->setGoal(yg);

  // simulate
  while (true)
  {
    // data logging
    Time = arma::join_horiz(Time, arma::vec({t}));
    Y_data = arma::join_horiz(Y_data, y);
    dY_data = arma::join_horiz(dY_data, dy);
    ddY_data = arma::join_horiz(ddY_data, ddy);
    // x_data = arma::join_horiz(x_data, arma::vec({x}));

    // gmp simulation
    arma::vec yc = arma::vec({0.0});
    arma::vec zc = arma::vec().zeros(Dim);
    gmp->update(s, y, z, yc, zc);
    dy = gmp->getYdot();
    dz = gmp->getZdot();
    arma::vec yc_dot = arma::vec({0.0});
    ddy = gmp->getYddot(yc_dot);

    // Stopping criteria
    if (t>=t_end) break;

    // Numerical integration
    t = t + dt;
    x = x + x_dot*dt;
    y = y + dy*dt;
    z = z + dz*dt;

    s.x = x;
  }

}

/*
void simulateGMPo_in_log_space(std::shared_ptr<gmp_::GMPo> &gmp_o, const arma::vec &Q0, const arma::vec &Qg,
                               double T, double dt, arma::rowvec &Time, arma::mat &Q_data, arma::mat &rotVel_data, arma::mat &rotAccel_data)
{
  // set initial values
  double t_end = T;
  double tau = t_end;

  double t = 0.0;
  double x = 0.0;
  double x_dot = 1/tau;
  double x_ddot = 0;
  arma::vec Q = Q0;
  arma::vec rotVel = arma::vec().zeros(3);
  arma::vec rotAccel = arma::vec().zeros(3);

  gmp_o->setQ0(Q0);  // set initial orientation
  gmp_o->setQg(Qg);  // set target orientation

  // simulate
  while (true)
  {
    // data logging
    Time = arma::join_horiz(Time, arma::vec({t}));
    Q_data = arma::join_horiz(Q_data, Q);
    rotVel_data = arma::join_horiz(rotVel_data, rotVel);
    rotAccel_data = arma::join_horiz(rotAccel_data, rotAccel);

    // GMP simulation
    arma::vec yc = arma::vec().zeros(3); // optional coupling for 'y' state
    arma::vec zc = arma::vec().zeros(3); // optional coupling for 'z' state
    arma::vec yc_dot = arma::vec().zeros(3); // derivative of coupling for 'y' state
    gmp_::Phase s(x,x_dot,x_ddot);
    rotAccel = gmp_o->calcRotAccel(s, Q, rotVel, Qg, yc, zc, yc_dot);

    // Stopping criteria
    if (t>1.5*t_end)
    {
      io_::PRINT_WARNING_MSG("Time limit reached... Stopping simulation!");
      break;
    }

    arma::vec eo = gmp_::quatLog(gmp_::quatProd(Qg, gmp_::quatInv(Q)));
    if (t>=t_end && arma::norm(eo)<0.02) break;

    // Numerical integration
    t = t + dt;
    x = x + x_dot*dt;
    Q = gmp_::quatProd( gmp_::quatExp(rotVel*dt), Q);
    rotVel = rotVel + rotAccel*dt;
  }

}


void simulateGMPo_in_quat_space(std::shared_ptr<gmp_::GMPo> &gmp_o, const arma::vec &Q0, const arma::vec &Qg,
                               double T, double dt, arma::rowvec &Time, arma::mat &Q_data, arma::mat &rotVel_data, arma::mat &rotAccel_data)
{
  // set initial values
  double t_end = T;
  double tau = t_end;

  double t = 0.0;
  double x = 0.0;
  double x_dot = 1/tau;
  double x_ddot = 0;
  arma::vec Q = Q0;
  arma::vec rotVel = arma::vec().zeros(3);
  arma::vec rotAccel = arma::vec().zeros(3);

  gmp_o->setQ0(Q0);  // set initial orientation
  gmp_o->setQg(Qg);  // set target orientation

  // simulate
  while (true)
  {
     // data logging
     Time = arma::join_horiz(Time, arma::vec({t}));
     Q_data = arma::join_horiz(Q_data, Q);
     rotVel_data = arma::join_horiz(rotVel_data, rotVel);
     rotAccel_data = arma::join_horiz(rotAccel_data, rotAccel);

     // GMP simulation
     arma::vec yc = arma::vec().zeros(3); // optional coupling for 'y' state
     arma::vec zc = arma::vec().zeros(3); // optional coupling for 'z' state
     arma::vec yc_dot = arma::vec().zeros(3); // derivative of coupling for 'y' state
     gmp_::Phase s(x,x_dot,x_ddot);
     rotAccel = gmp_o->calcRotAccel(s, Q, rotVel, Qg, yc, zc, yc_dot);

     // Stopping criteria
     if (t>1.5*t_end)
     {
       io_::PRINT_WARNING_MSG("Time limit reached... Stopping simulation!");
       break;
     }

     arma::vec eo = gmp_::quatLog(gmp_::quatProd(Qg, gmp_::quatInv(Q)));
     if (t>=t_end && arma::norm(eo)<0.02) break;

     // Numerical integration
     t = t + dt;
     x = x + x_dot*dt;
     Q = gmp_::quatProd( gmp_::quatExp(rotVel*dt), Q);
     rotVel = rotVel + rotAccel*dt;

  }

}

*/
