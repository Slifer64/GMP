#include <gmp_lib/GMP/GMP_nDoF_IO.h>

namespace as64_
{

namespace gmp_
{

  GMP_nDoF_IO::GMP_nDoF_IO(const std::shared_ptr<gmp_::GMP_nDoF> &gmp)
  {
    this->gmp = gmp;
  }


  void GMP_nDoF_IO::write(const std::string &filename, const std::string &prefix)
  {
    gmp_::FileIO fid(filename, gmp_::FileIO::out|gmp_::FileIO::trunc);
    write(fid, prefix);
  }

  void GMP_nDoF_IO::write(gmp_::FileIO &fid, const std::string &prefix)
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

  void GMP_nDoF_IO::read(const std::string &filename, const std::string &prefix)
  {
    gmp_::FileIO fid(filename, gmp_::FileIO::in);
    read(fid, prefix);
  }

  void GMP_nDoF_IO::read(gmp_::FileIO &fid, const std::string &prefix)
  {
    fid.read(prefix + "weights", gmp->W);
    fid.read(prefix + "damping", gmp->D);
    fid.read(prefix + "stiffness", gmp->K);
    //fid.read(prefix + "N_kernels", N_kernels);
    //fid.read(prefix + "N_DoFs", n_dofs);
    int scale_type;
    fid.read(prefix + "scale_type", scale_type);
    fid.read(prefix + "c", gmp->c);
    fid.read(prefix + "h", gmp->h);

    gmp->Y0d = gmp->W*gmp->regressVec(0);
    gmp->Ygd = gmp->W*gmp->regressVec(1);

    gmp->setY0(gmp->Y0d);
    gmp->setGoal(gmp->Ygd);

    unsigned n_dofs = gmp->numOfDoFs();
    gmp->y_dot = arma::vec().zeros(n_dofs);
    gmp->z_dot = arma::vec().zeros(n_dofs);

    if (scale_type == TrajScale::PROP_SCALE) gmp->setScaleMethod( gmp_::TrajScale::Ptr(new gmp_::TrajScale_Prop(n_dofs)) );
    else if (scale_type == TrajScale::ROT_MIN_SCALE) gmp->setScaleMethod( gmp_::TrajScale::Ptr(new gmp_::TrajScale_Rot_min()) );
    else if (scale_type == TrajScale::ROT_WB_SCALE) gmp->setScaleMethod( gmp_::TrajScale::Ptr(new gmp_::TrajScale_Rot_wb()) );
    else throw std::runtime_error("[GMP_nDoF_IO::read]: Unsupported scale type \"" + std::to_string(scale_type) + "\"...\n");
  }


} // namespace gmp_

} // namespace as64_
