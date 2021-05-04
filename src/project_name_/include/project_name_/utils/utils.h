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

void throw_error(const std::string &err_msg);


void simulateGMP_nDoF(gmp_::GMP_nDoF::Ptr &gmp, const arma::vec &y0,
                      const arma::vec &yg, double T, double dt,
                      arma::mat &Time, arma::mat &Y_data, arma::mat &dY_data, arma::mat &ddY_data);

void simulateGMPo_in_log_space(gmp_::GMPo::Ptr &gmp_o, const arma::vec &Q0, const arma::vec &Qg,
                               double T, double dt, arma::rowvec &Time, arma::mat &Q_data, arma::mat &rotVel_data, arma::mat &rotAccel_data);

void simulateGMPo_in_quat_space(gmp_::GMPo::Ptr &gmp_o, const arma::vec &Q0, const arma::vec &Qg,
                               double T, double dt, arma::rowvec &Time, arma::mat &Q_data, arma::mat &rotVel_data, arma::mat &rotAccel_data);

void simulateGMPo_in_Cart_space(gmp_::GMPo::Ptr &gmp_o, const arma::vec &Q0, const arma::vec &Qg,
                               double T, double dt, arma::rowvec &Time, arma::mat &Q_data, arma::mat &rotVel_data, arma::mat &rotAccel_data);

#endif // GMP_TEST_UTILS_H
