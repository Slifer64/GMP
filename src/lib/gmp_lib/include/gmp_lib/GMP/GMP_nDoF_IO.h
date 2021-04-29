#ifndef GMP_LIB_GMP_NDOF_IO_H
#define GMP_LIB_GMP_NDOF_IO_H

#include <gmp_lib/GMP/GMP_nDoF.h>
#include <gmp_lib/io/file_io.h>

namespace as64_
{

namespace gmp_
{

class GMP_nDoF_IO
{

public:

  typedef std::shared_ptr<GMP_nDoF_IO> Ptr;

  /** Constructor.
   *  @param[in] gmp: n_DoF dmp.
   */
  GMP_nDoF_IO(const std::shared_ptr<gmp_::GMP_nDoF> &gmp);

  /** Write the GMP model to a file.
   * @param[in] filename: String with the name of the file.
   * @param[in] prefix: The prefix that will be used for writing the names of all GMP params (optional, default="").
   */
  void write(const std::string &filename, const std::string &prefix="");

  /** Write the GMP model to a file.
   * @param[in] fid: Object of type @FileIO associated with the file.
   * @param[in] prefix: The prefix that will be used for writing the names of all GMP params (optional, default="").
   */
  void write(gmp_::FileIO &fid, const std::string &prefix="");

  /** Reads the GMP model from a file.
   * @param[in] fid: Filename string or object of type @FileIO associated with the file.
   * @param[in] prefix: The prefix that will be used for reading the names of all GMP params (optional, default="").
   */
  void read(const std::string &filename, const std::string &prefix="");

  /** Reads the GMP model from a file.
   * @param[in] fid: Filename string or object of type @FileIO associated with the file.
   * @param[in] prefix: The prefix that will be used for reading the names of all GMP params (optional, default="").
   */
  void read(gmp_::FileIO &fid, const std::string &prefix="");

private:

  std::shared_ptr<gmp_::GMP_nDoF> gmp;

}; // GMP_nDoF_IO

} // namespace gmp_

} // namespace as64_

#endif // GMP_LIB_GMP_NDOF_IO_H
