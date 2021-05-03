#include <gmp_lib/GMP/GMPo_IO.h>
#include <gmp_lib/GMP/GMP_nDoF_IO.h>

namespace as64_
{

namespace gmp_
{

  void GMPo_IO::write(const gmp_::GMPo *gmp, const std::string &filename, const std::string &prefix)
  {
    gmp_::FileIO fid(filename, gmp_::FileIO::out|gmp_::FileIO::trunc);
    write(gmp, fid, prefix);
  }

  void GMPo_IO::write(const gmp_::GMPo *gmp, gmp_::FileIO &fid, const std::string &prefix)
  {
    GMP_nDoF_IO::write(dynamic_cast<const GMP_nDoF *>(gmp), fid, prefix);
    fid.write(prefix + "Qd0", gmp->Qd0);
  }

  void GMPo_IO::read(gmp_::GMPo *gmp, const std::string &filename, const std::string &prefix)
  {
    gmp_::FileIO fid(filename, gmp_::FileIO::in);
    read(gmp, fid, prefix);
  }

  void GMPo_IO::read(gmp_::GMPo *gmp, gmp_::FileIO &fid, const std::string &prefix)
  {
    GMP_nDoF_IO::read(dynamic_cast<GMP_nDoF *>(gmp), fid, prefix);
    fid.read(prefix + "Qd0", gmp->Qd0);
    gmp->setQ0(gmp->Qd0);
  }


} // namespace gmp_

} // namespace as64_
