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
// #include <gmp_lib/math/quaternions.h>

using namespace as64_;

void PRINT_INFO_MSG(const std::string &msg, std::ostream &out = std::cout);
void PRINT_CONFIRM_MSG(const std::string &msg, std::ostream &out = std::cout);
void PRINT_WARNING_MSG(const std::string &msg, std::ostream &out = std::cout);
void PRINT_ERROR_MSG(const std::string &msg, std::ostream &out = std::cout);

void simulateGMP_nDoF(gmp_::GMP_nDoF::Ptr &gmp, const arma::vec &y0,
                      const arma::vec &yg, double T, double dt,
                      arma::mat &Time, arma::mat &Y_data, arma::mat &dY_data, arma::mat &ddY_data);

/*
void simulateGMPo_in_log_space(std::shared_ptr<gmp_::GMPo> &gmp_o, const arma::vec &Q0, const arma::vec &Qg,
                               double T, double dt, arma::rowvec &Time, arma::mat &Q_data, arma::mat &rotVel_data, arma::mat &rotAccel_data);


void simulateGMPo_in_quat_space(std::shared_ptr<gmp_::GMPo> &gmp_o, const arma::vec &Q0, const arma::vec &Qg,
                              double T, double dt, arma::rowvec &Time, arma::mat &Q_data, arma::mat &rotVel_data, arma::mat &rotAccel_data);
*/

#endif // GMP_TEST_UTILS_H
