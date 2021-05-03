#ifndef GMP_TEST_UTILS_H
#define GMP_TEST_UTILS_H

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>
#include <memory>
#include <armadillo>
#include <gmp_lib/gmp_lib.h>
#include <gmp_lib/math/quaternions.h>

#define PROJECT_NAME_ "project_name_"

using namespace as64_;

class Timer
{
public:
  static void tic() { Timer::timer.tic(); }
  static void toc() { std::cerr << "Elapsed time " + std::to_string(Timer::timer.toc()) + " sec\n"; }
private:
  static arma::wall_clock timer;
};

void PRINT_INFO_MSG(const std::string &msg, std::ostream &out = std::cout);
void PRINT_CONFIRM_MSG(const std::string &msg, std::ostream &out = std::cout);
void PRINT_WARNING_MSG(const std::string &msg, std::ostream &out = std::cout);
void PRINT_ERROR_MSG(const std::string &msg, std::ostream &out = std::cout);

void simulateGMP_nDoF(gmp_::GMP_nDoF::Ptr &gmp, const arma::vec &y0,
                      const arma::vec &yg, double T, double dt,
                      arma::mat &Time, arma::mat &Y_data, arma::mat &dY_data, arma::mat &ddY_data);
/*
class SimOptions
{
public:

  enum StopCriteria
  {
    TIME_DURATION = 1,
    POS_ERROR_THRES = 1<<1,
  }

  enum ExitStatus
  {
    TIME_LIMIT_EXCEEDED = 0,
    STOP_CRITERIA_SATIFIED = 1
  }

  SimOptions(double time_duration, double pos_err_threshold, double time_limit, int stop_criteria=TIME_DURATION|POS_ERROR_THRES)
  {
    this->time_duration = time_duration;
    this->pos_err_threshold = pos_err_threshold;
    this->time_limit = time_limit;
    this->stop_criteria = stop_criteria;
  }

  bool finish(double pos_err, double t) const
  {
    bool success = true;
    if (stop_criteria&TIME_DURATION) success &= (t<time_duration);
    if (stop_criteria&POS_ERROR_THRES) success &= (pos_err<pos_err_threshold);

    if (success)
    {
      exit_status = STOP_CRITERIA_SATIFIED;
      return true;
    }
    else if (t > time_limit)
    {
      exit_status = TIME_LIMIT_EXCEEDED;
      return true;
    }
    else return false;
  }

  ExitStatus exit_status;

private:

  double time_limit;
  double time_duration;
  double sim_timestep;
  double pos_err_threshold;

  int stop_criteria;
};*/

void simulateGMPo_in_log_space(gmp_::GMPo::Ptr &gmp_o, const arma::vec &Q0, const arma::vec &Qg,
                               double T, double dt, arma::rowvec &Time, arma::mat &Q_data, arma::mat &rotVel_data, arma::mat &rotAccel_data);

void simulateGMPo_in_quat_space(gmp_::GMPo::Ptr &gmp_o, const arma::vec &Q0, const arma::vec &Qg,
                               double T, double dt, arma::rowvec &Time, arma::mat &Q_data, arma::mat &rotVel_data, arma::mat &rotAccel_data);

void simulateGMPo_in_Cart_space(gmp_::GMPo::Ptr &gmp_o, const arma::vec &Q0, const arma::vec &Qg,
                               double T, double dt, arma::rowvec &Time, arma::mat &Q_data, arma::mat &rotVel_data, arma::mat &rotAccel_data);

#endif // GMP_TEST_UTILS_H
