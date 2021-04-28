#ifndef AS64_gmp_IO_H
#define AS64_gmp_IO_H

#include <iostream>
#include <fstream>
#include <memory>
#include <iomanip>
#include <cerrno>
#include <streambuf>
#include <armadillo>

namespace as64_
{

namespace gmp_
{

/** \brief Reads a 2D matrix from stream \a in.
 *  \details Reads the elements of the matrix row by row.
 *  @param[out] m The 2D matrix
 *  @param[in] n_rows Number of rows of the matrix
 *  @param[in] n_cols Number of columns of the matrix
 *  @param[in] in The input stream (default = std::cin).
 *  @param[in] binary Flag indicating the format (true for binary, false for text, optional, default = false).
 */
void read_mat(arma::mat &m, long n_rows, long n_cols, std::istream &in = std::cin, bool binary = true);


/** \brief Reads a 2D matrix from stream \a in.
 *  \details Reads first the number of rows and columns and then the elements of the matrix row by row.
 *  @param[out] m The 2D matrix
 *  @param[in] in The input stream (default = std::cin).
 *  @param[in] binary Flag indicating the format (true for binary, false for text, optional, default = false).
 */
void read_mat(arma::mat &m, std::istream &in = std::cin, bool binary = true);


/** \brief Reads a scalar value from stream \a in.
 *  \details Reads a scalar value in the format specified by \a binary flag.
 *  @param[out] scalar scalar value.
 *  @param[in] in The input stream (default = std::cin).
 *  @param[in] binary Flag indicating the format (true for binary, false for text, optional, default = false).
 */
template <typename T>
void read_scalar(T &scalar, std::istream &in = std::cin, bool binary = true)
{
  if (binary) in.read((char *)(&scalar), sizeof(scalar));
  else in >> scalar;
}

/** \brief Writes a 2D matrix in stream 'fid'.
 *  \details Writes the elements of the matrix row by row.
 *  @param[in] m The 2D matrix
 *  @param[in] n_rows Number of rows of the matrix
 *  @param[in] n_cols Number of columns of the matrix
 *  @param[in] out The output stream (optional, default = 1 for output to screen).
 *  @param[in] binary Flag indicating the format (true for binary, false for text, optional, default = false).
 *  @param[in] precision precision in txt format (optional, default = 6).
 */
void write_mat(const arma::mat &m, int n_rows, int n_cols, std::ostream &out = std::cout, bool binary = true, int precision = 6);


/** \brief Writes a 2D matrix in stream \a out.
 *  \details Writes first the number of rows and columns and then the elements of the matrix row by row.
 *  @param[in] m The 2D matrix
 *  @param[in] out The output stream (optional, default = 1 for output to screen).
 *  @param[in] binary Flag indicating the format (true for binary, false for text, optional, default = false).
 *  @param[in] precision precision in txt format (optional, default = 6).
 */
void write_mat(const arma::mat &m, std::ostream &out = std::cout, bool binary = true, int precision = 6);


/** \brief Writes a scalar value in stream \a out
 *  \details Writes a scalar in the format specified by \a binary flag. If the format is binary,
 *           the class type of \a scalar is used to determine the number of bytes to use.
 *  @param[in] scalar scalar value.
 *  @param[in] out The output stream (optional, default = std::cout).
 *  @param[in] binary Flag indicating the format (true for binary, false for text, optional, default = false).
 *  @param[in] precision precision in txt format (optional, default = 6).
 */
template <typename T>
void write_scalar(T scalar, std::ostream &out = std::cout, bool binary = true, int precision = 6)
{
  if (binary) out.write((const char *)(&scalar), sizeof(scalar));
  else out << std::setprecision(precision) << scalar << "\n";
}

} // namespace gmp_

} // namespace as64_


#endif // AS64_gmp_IO_H
