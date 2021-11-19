#include <gmp_test/utils/utils.h>

arma::wall_clock Timer::timer;


void PRINT_INFO_MSG(const std::string &msg, std::ostream &out)
{
  out << "\033[1m" << "\033[34m" << msg << "\033[0m";
}

void PRINT_CONFIRM_MSG(const std::string &msg, std::ostream &out)
{
  out << "\033[1m" << "\033[32m" << msg << "\033[0m";
}

void PRINT_WARNING_MSG(const std::string &msg, std::ostream &out)
{
  out << "\033[1m" << "\033[33m" << msg << "\033[0m";
}

void PRINT_ERROR_MSG(const std::string &msg, std::ostream &out)
{
  out << "\033[1m" << "\033[31m" << msg << "\033[0m";
}

void throw_error(const std::string &err_msg)
{
  PRINT_ERROR_MSG(err_msg);
  exit(-1);
}

void simulateGMP(gmp_::GMP::Ptr &gmp, const arma::vec &y0,
                 const arma::vec &yg, double T, double dt,
                 arma::mat &Time, arma::mat &Y_data, arma::mat &dY_data, arma::mat &ddY_data)
{

  int Dim = gmp->numOfDoFs();
  arma::vec y = y0;
  arma::vec dy = arma::vec().zeros(Dim);
  arma::vec ddy = arma::vec().zeros(Dim);
  // arma::vec z = arma::vec().zeros(Dim);
  // arma::vec dz = arma::vec().zeros(Dim);

  double t = 0.0;
  double t_end = T;
  double tau = t_end;
  double x = 0;
  double x_dot = 1/tau;
  double x_ddot = 0;
  // gmp_::Phase s(x, x_dot, x_ddot);

  gmp->setY0(y0); // set initial position
  gmp->setGoal(yg); // set target/final position

  // simulate
  while (true)
  {
    // data logging
    Time = arma::join_horiz(Time, arma::vec({t}));
    Y_data = arma::join_horiz(Y_data, y);
    dY_data = arma::join_horiz(dY_data, dy);
    ddY_data = arma::join_horiz(ddY_data, ddy);
    // x_data = arma::join_horiz(x_data, arma::vec({x}));

    // Update the target if you want...
    // gmp->setGoal(g);

    // gmp simulation

    // Deprecated, old DMP style evolution:
    // arma::vec yc = arma::vec({0.0});
    // arma::vec zc = arma::vec().zeros(Dim);
    // gmp->update(s, y, z, yc, zc);
    // dy = gmp->getYdot();
    // dz = gmp->getZdot();
    // arma::vec yc_dot = arma::vec({0.0});
    // ddy = gmp->getYddot(yc_dot);

    // get the scaled DMP learned trajectory
    arma::vec y_x = gmp->getYd(x);
    arma::vec dy_x = gmp->getYdDot(x, x_dot);
    arma::vec ddy_x = gmp->getYdDDot(x, x_dot, x_ddot);
    
    double K = 300; // set the DMP stiffness
    double D = 60; // set the DMP damping
    
    arma::vec external_signal = arma::vec().zeros(Dim); // optionally add some external signal
    
    // Track it using a 2nd order dynamical system. This is actually the DMP. 
    ddy = ddy_x + D*(dy_x - dy) + K*(y_x - y) + external_signal;

    // if external_signal == 0 you can obviously set directly:
//    y = y_x;
//    dy = dy_x;
//    ddy = ddy_x;

    // Stopping criteria
    double ep = arma::norm(y-yg);
    if (t>=t_end && ep<1e-3) // && arma::norm(dy)<5e-3)
    {
      PRINT_CONFIRM_MSG("Target reached!\nDistance from target: " + std::to_string(ep) + "\n");
      break;
    }

    if (t>=1.05*t_end)
    {
      PRINT_WARNING_MSG("Time limit exceeded...\nDistance from target: " + std::to_string(ep) + "\nStopping simulation!\n");
      break;
    }

    // Numerical integration
    t = t + dt;
    x = x + x_dot*dt;
    x_dot = x_dot + x_ddot*dt;
    y = y + dy*dt;
    dy = dy + ddy*dt;
    //z = z + dz*dt;

    //s.x = x;
  }

}


/** Simulates a dmp encoding Cartesian orientation using unit quaternions.
 * The DMP dynamics are formulated in the Cartesian orientation space.
 * All quaternions are assumed to be unit, with format [w; x; y; z]
 * @param[in] gmp: object of type @GMPo.
 * @param[in] Q0: Initial orientation as a unit quaternion.
 * @param[in] Qg: Target orientation as a unit quaternion.
 * @param[in] T: Movement total time duration.
 * @param[in] dt: Simulation timestep.
 * @param[out] Time: 1 x N rowvector with simulation output timestamps.
 * @param[out] Q_data: 4 x N matrix with simulation output positions.
 * @param[out] rotVel_data: 3 x N matrix with simulation output rotational velocities.
 * @param[out] rotAccel_data: 3 x N matrix with simulation output rotational accelerations.
 */
void simulateGMPo_in_Cart_space(gmp_::GMPo::Ptr &gmp_o, const arma::vec &Q0, const arma::vec &Qg,
                               double T, double dt, arma::rowvec &Time, arma::mat &Q_data, arma::mat &rotVel_data, arma::mat &rotAccel_data)
{
  // set initial values
  double t_end = T;
  double tau = t_end;

  double t = 0.0;

  // phase variable, from 0 to 1
  double x = 0.0;
  double x_dot = 1/tau;
  double x_ddot = 0;

  // Cartesian orientation trajectory
  arma::vec Q = Q0;
  arma::vec rotVel = arma::vec().zeros(3);
  arma::vec rotAccel = arma::vec().zeros(3);

  gmp_o->setQ0(Q0);   // set initial orientation
  gmp_o->setQg(Qg);   // set target orientation

  // simulate
  while (true)
  {
      // data logging
      Time = arma::join_horiz(Time, arma::vec({t}));
      Q_data = arma::join_horiz(Q_data, Q);
      rotVel_data = arma::join_horiz(rotVel_data, rotVel);
      rotAccel_data = arma::join_horiz(rotAccel_data, rotAccel);

      // GMP simulation
      arma::vec Qd = gmp_o->getQd(x);
      arma::vec Vd = gmp_o->getVd(x, x_dot);
      arma::vec Vd_dot = gmp_o->getVdDot(x, x_dot, x_ddot);
      // Or get them all together:
      // gmp_o->getRefTraj(x, x_dot, x_ddot, Qd, Vd, Vd_dot);

      arma::vec external_signal = arma::vec().zeros(3); // optionally add some external signal

      rotAccel = Vd_dot + 5*(Vd-rotVel) + 20*gmp_::quatLog(gmp_::quatDiff(Qd,Q)) + external_signal;

      // if external_signal == 0 you can obviously set directly:
//     Q = Qd;
//     rotVel = Vd;
//     rotAccel = Vd_dot;

      // Stopping criteria
      arma::vec eo = gmp_::quatLog(gmp_::quatProd(Qg, gmp_::quatInv(Q)));
      if (t>=t_end && arma::norm(eo)<0.01)
      {
        PRINT_CONFIRM_MSG("Target reached!\nDistance from target: " + std::to_string(arma::norm(eo)) + "\n");
        break;
      }

      if (t>=1.05*t_end)
      {
        PRINT_WARNING_MSG("Time limit exceeded...\nDistance from target: eo = " + std::to_string(arma::norm(eo)) + "\nStopping simulation!\n");
        break;
      }

      // Numerical integration
      t = t + dt;
      x = x + x_dot*dt;
      x_dot = x_dot + x_ddot*dt;
      Q = gmp_::quatProd( gmp_::quatExp(rotVel*dt), Q);
      rotVel = rotVel + rotAccel*dt;
  }

}

/** Simulates a dmp encoding Cartesian orientation using unit quaternions.
 * The DMP dynamics are formulated in the quaterion-logarith space and the
 * generated trajectory is converted to the Cartesian orientation space.
 * All quaternions are assumed to be unit, with format [w; x; y; z]
 * @param[in] gmp: object of type @GMPo.
 * @param[in] Q0: Initial orientation as a unit quaternion.
 * @param[in] Qg: Target orientation as a unit quaternion.
 * @param[in] T: Movement total time duration.
 * @param[in] dt: Simulation timestep.
 * @param[out] Time: 1 x N rowvector with simulation output timestamps.
 * @param[out] Q_data: 4 x N matrix with simulation output positions.
 * @param[out] rotVel_data: 3 x N matrix with simulation output rotational velocities.
 * @param[out] rotAccel_data: 3 x N matrix with simulation output rotational accelerations.
 */
void simulateGMPo_in_log_space(gmp_::GMPo::Ptr &gmp_o, const arma::vec &Q0, const arma::vec &Qg,
                               double T, double dt, arma::rowvec &Time, arma::mat &Q_data, arma::mat &rotVel_data, arma::mat &rotAccel_data)
{
  // set initial values
  double t_end = T;
  double tau = t_end;

  double t = 0.0;

  // phase variable, from 0 to 1
  double x = 0.0;
  double x_dot = 1/tau;
  double x_ddot = 0;

  // Cartesian orientation trajectory
  arma::vec Q = Q0;
  arma::vec Q_prev = Q;
  arma::vec rotVel = arma::vec().zeros(3);
  arma::vec rotAccel = arma::vec().zeros(3);

  // quaternion logarithm trajectory
  arma::vec y = gmp_::GMPo::quat2q(Q0, Q0);
  arma::vec dy = arma::vec().zeros(3);
  arma::vec ddy = arma::vec().zeros(3);

  gmp_o->setQ0(Q0); // set initial orientation
  gmp_o->setQg(Qg); // set target orientation

  // simulate
  while (true)
  {
    // data logging
    Time = arma::join_horiz(Time, arma::vec({t}));
    Q_data = arma::join_horiz(Q_data, Q);
    rotVel_data = arma::join_horiz(rotVel_data, rotVel);
    rotAccel_data = arma::join_horiz(rotAccel_data, rotAccel);

    // DMP simulation
    
    // get the scaled DMP learned trajectory in the quatLog space.
    arma::vec y_x = gmp_o->getYd(x);
    arma::vec dy_x = gmp_o->getYdDot(x, x_dot);
    arma::vec ddy_x = gmp_o->getYdDDot(x, x_dot, x_ddot);
    
    double K = 20; // set the DMP stiffness in the quatLog space
    double D = 5; // set the DMP damping in the quatLog space
    
    arma::vec external_signal = arma::vec().zeros(3); // optionally add some external signal
    
    // Track it using a 2nd order dynamical system.
    // DMP dynamics in the quatLog space:
    ddy = ddy_x + D*(dy_x - dy) + K*(y_x - y) + external_signal;
    
    // if external_signal == 0 you can obviously set directly:
//     y = y_x;
//     dy = dy_x;
//     ddy = ddy_x;

    // Stopping criteria
    arma::vec eo = gmp_::quatLog(gmp_::quatProd(Qg, gmp_::quatInv(Q)));
    if (t>=t_end && arma::norm(eo)<0.01)
    {
      PRINT_CONFIRM_MSG("Target reached!\nDistance from target: " + std::to_string(arma::norm(eo)) + "\n");
      break;
    }

    if (t>=1.05*t_end)
    {
      PRINT_WARNING_MSG("Time limit exceeded...\nDistance from target: eo = " + std::to_string(arma::norm(eo)) + "\nStopping simulation!\n");
      break;
    }

    // Numerical integration
    t = t + dt;
    x = x + x_dot*dt;
    x_dot = x_dot + x_ddot*dt;
    y = y + dy*dt;
    dy = dy + ddy*dt;

  
    // Convert trajectory from quatLog to the Cartesian orientation space
    Q = gmp_o->q2quat(y, Q0);
    // This is required to avoid discontinuities (since Q and -Q are the same orientation)
    if (arma::dot(Q_prev,Q)<0) Q = -Q; 
    Q_prev = Q;
    
    arma::vec Q1 = gmp_o->getQ1(Q, Q0);
    rotVel = gmp_::qLogDot_to_rotVel(dy, Q1);
    rotAccel = gmp_::qLogDDot_to_rotAccel(ddy, rotVel, Q1);

  }

}
