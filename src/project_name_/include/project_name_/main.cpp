#include <iostream>
#include <memory>
#include <ros/ros.h>

#include <gmp_lib/gmp_lib.h>

using namespace as64_;

int main(int argc, char **argv)
{
  ros::init(argc, argv, "$project_name$");

  gmp_::GMP_nDoF::Ptr gmp(new gmp_::GMP_nDoF(1,2));

  gmp_::GMP_nDoF_IO::Ptr gmp_io(new gmp_::GMP_nDoF_IO(gmp));

  gmp_::GMP_nDoF_Update::Ptr gmp_up(new gmp_::GMP_nDoF_Update(gmp));


  return 0;
}
