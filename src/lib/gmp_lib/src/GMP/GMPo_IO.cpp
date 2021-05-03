#include <gmp_lib/GMP/GMPo_IO.h>
#include <gmp_lib/GMP/GMP_nDoF_IO.h>

namespace as64_
{

namespace gmp_
{

  void GMPo_IO::write(const gmp_::GMPo::Ptr gmp, const std::string &filename, const std::string &prefix)
  {
    gmp_::FileIO fid(filename, gmp_::FileIO::out|gmp_::FileIO::trunc);
    write(gmp, fid, prefix);
  }

  void GMPo_IO::write(const gmp_::GMPo::Ptr gmp, gmp_::FileIO &fid, const std::string &prefix)
  {
    /*
    unsigned N_kernels = gmp->numOfKernels();
    unsigned n_dofs = gmp->numOfDoFs();

    fid.write(prefix + "weights", gmp->W);
    fid.write(prefix + "damping", gmp->D);
    fid.write(prefix + "stiffness", gmp->K);
    fid.write(prefix + "N_kernels", N_kernels);
    fid.write(prefix + "N_DoFs", n_dofs);
    fid.write(prefix + "scale_type", (int)(gmp->traj_sc->getScaleType()));
    fid.write(prefix + "c", gmp->c);
    fid.write(prefix + "h", gmp->h);
    */
  }

  void GMPo_IO::read(gmp_::GMPo::Ptr gmp, const std::string &filename, const std::string &prefix)
  {
    gmp_::FileIO fid(filename, gmp_::FileIO::in);
    read(gmp, fid, prefix);
  }

  void GMPo_IO::read(gmp_::GMPo::Ptr gmp, gmp_::FileIO &fid, const std::string &prefix)
  {
    /*
    fid.read(prefix + "weights", gmp->W);
    fid.read(prefix + "damping", gmp->D);
    fid.read(prefix + "stiffness", gmp->K);
    //fid.read(prefix + "N_kernels", N_kernels);
    //fid.read(prefix + "N_DoFs", n_dofs);
    int scale_type;
    fid.read(prefix + "scale_type", scale_type);
    fid.read(prefix + "c", gmp->c);
    fid.read(prefix + "h", gmp->h);

    // std::cerr << "w: " << gmp->W.n_rows << " x  " << gmp->W.n_cols << "\n";
    // std::cerr << "D: " << gmp->D.n_rows << " x  " << gmp->D.n_cols << "\n";
    // std::cerr << "K: " << gmp->K.n_rows << " x  " << gmp->K.n_cols << "\n";
    // std::cerr << "c: " << gmp->c.n_rows << " x  " << gmp->c.n_cols << "\n";
    // std::cerr << "h: " << gmp->h.n_rows << " x  " << gmp->h.n_cols << "\n";

    unsigned n_dofs = gmp->numOfDoFs();

    gmp->Y0d = gmp->W*gmp->regressVec(0);
    gmp->Ygd = gmp->W*gmp->regressVec(1);
    gmp->Y0 = gmp->Y0d;
    gmp->Yg = gmp->Ygd;
    gmp->y_dot = arma::vec().zeros(n_dofs);
    gmp->z_dot = arma::vec().zeros(n_dofs);

    if (scale_type == TrajScale::PROP_SCALE) gmp->setScaleMethod( gmp_::TrajScale::Ptr(new gmp_::TrajScale_Prop(n_dofs)) );
    else if (scale_type == TrajScale::ROT_MIN_SCALE) gmp->setScaleMethod( gmp_::TrajScale::Ptr(new gmp_::TrajScale_Rot_min()) );
    else if (scale_type == TrajScale::ROT_WB_SCALE) gmp->setScaleMethod( gmp_::TrajScale::Ptr(new gmp_::TrajScale_Rot_wb()) );
    else throw std::runtime_error("[GMPo_IO::read]: Unsupported scale type \"" + std::to_string(scale_type) + "\"...\n");
    */
  }


} // namespace gmp_

} // namespace as64_
