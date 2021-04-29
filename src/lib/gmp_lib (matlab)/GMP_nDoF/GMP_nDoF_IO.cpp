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

  /** Constructor.
   *  @param[in] gmp: n_DoF dmp.
   */
  GMP_nDoF_IO(const std::shared_ptr<gmp_::GMP_nDoF> &gmp)
  {
    this->gmp = gmp;
  }

  /** Write the GMP model to a file.
   * @param[in] filename: String with the name of the file.
   * @param[in] prefix: The prefix that will be used for writing the names of all GMP params (optional, default="").
   */
  void write(const std::string &filename, const std::string &prefix="")
  {
    gmp_::FileIO fid(filename, gmp_::FileIO::out|gmp_::FileIO::trunc);
    write(fid, prefix);
  }

  /** Write the GMP model to a file.
   * @param[in] fid: Object of type @FileIO associated with the file.
   * @param[in] prefix: The prefix that will be used for writing the names of all GMP params (optional, default="").
   */
  void write(gmp_::FileIO &fid, const std::string &prefix="")
  {
    unsigned N_kernels = gmp->numOfKernels();
    unsigned n_dofs = gmp->numOfDoFs();

    fid.write(prefix + "weights", gmp->W);
    fid.write(prefix + "damping", gmp->D);
    fid.write(prefix + "stiffness", gmp->K);
    fid.write(prefix + "N_kernels", N_kernels);
    fid.write(prefix + "N_DoFs", n_dofs);
    fid.write(prefix + "scale_type", gmp->traj_sc->getScaleType());
    fid.write(prefix + "c", gmp->c);
    fid.write(prefix + "h", gmp->h);
  }

  /** Reads the GMP model from a file.
   * @param[in] fid: Filename string or object of type @FileIO associated with the file.
   * @param[in] prefix: The prefix that will be used for reading the names of all GMP params (optional, default="").
   */
  void read(const std::string &filename, const std::string &prefix="")
  {
    gmp_::FileIOfid fid(filename, gmp_::FileIO::in);
    read(fid, prefix);
  }

  /** Reads the GMP model from a file.
   * @param[in] fid: Filename string or object of type @FileIO associated with the file.
   * @param[in] prefix: The prefix that will be used for reading the names of all GMP params (optional, default="").
   */
  void read(gmp_::FileIO &fid, const std::string &prefix="")
  {
    fid.read(prefix + "weights", gmp->W);
    fid.read(prefix + "damping", gmp->D);
    fid.read(prefix + "stiffness", gmp->K);
    //fid.read(prefix + "N_kernels", N_kernels);
    //fid.read(prefix + "N_DoFs", n_dofs);
    fid.read(prefix + "c", gmp->c);
    fid.read(prefix + "h", gmp->h);

    gmp->Y0d = gmp->W*gmp->regressVec(0);
    gmp->Ygd = gmp->W*gmp->regressVec(1);

    gmp->setY0(gmp->Y0d);
    gmp->setGoal(gmp->Ygd);

    unsigned n_dofs = gmp->numOfDoFs();
    gmp->y_dot = arma::vec().zeros(n_dofs);
    gmp->z_dot = arma::vec().zeros(n_dofs);
  }

private:

  std::shared_ptr<gmp_::GMP_nDoF> gmp;

}; // GMP_nDoF_IO

} // namespace gmp_

} // namespace as64_

#endif // GMP_LIB_GMP_NDOF_IO_H
