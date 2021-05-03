#ifndef GMP_LIB_GMP_ORIENT_IO_H
#define GMP_LIB_GMP_ORIENT_IO_H

#include <gmp_lib/GMP/GMPo.h>
#include <gmp_lib/io/file_io.h>

namespace as64_
{

namespace gmp_
{

class GMPo_IO
{

public:

  //GMPo_IO() {}

  /** Write the GMP model to a file.
   * @param[in] filename: String with the name of the file.
   * @param[in] prefix: The prefix that will be used for writing the names of all GMP params (optional, default="").
   */
  static void write(const gmp_::GMPo::Ptr gmp, const std::string &filename, const std::string &prefix="");

  /** Write the GMP model to a file.
   * @param[in] fid: Object of type @FileIO associated with the file.
   * @param[in] prefix: The prefix that will be used for writing the names of all GMP params (optional, default="").
   */
  static void write(const gmp_::GMPo::Ptr gmp, gmp_::FileIO &fid, const std::string &prefix="");

  /** Reads the GMP model from a file.
   * @param[in] fid: Filename string or object of type @FileIO associated with the file.
   * @param[in] prefix: The prefix that will be used for reading the names of all GMP params (optional, default="").
   */
  static void read(gmp_::GMPo::Ptr gmp, const std::string &filename, const std::string &prefix="");

  /** Reads the GMP model from a file.
   * @param[in] fid: Filename string or object of type @FileIO associated with the file.
   * @param[in] prefix: The prefix that will be used for reading the names of all GMP params (optional, default="").
   */
  static void read(gmp_::GMPo::Ptr gmp, gmp_::FileIO &fid, const std::string &prefix="");

}; // GMPo_IO

} // namespace gmp_

} // namespace as64_

#endif // GMP_LIB_GMP_ORIENT_IO_H
