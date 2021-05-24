#include <iostream>
#include <memory>
#include <ros/ros.h>

#include <gmp_lib/gmp_lib.h>

using namespace as64_;

int main(int argc, char **argv)
{
  ros::init(argc, argv, "$project_name$");

  gmp_::GMP::Ptr gmp(new gmp_::GMP(1,2));

  // gmp_::GMP_IO::Ptr gmp_io(new gmp_::GMP_IO(gmp));
  //
  // gmp_::GMP_Update::Ptr gmp_up(new gmp_::GMP_Update(gmp));

  std::cerr << "==========   Everything is ok!  ==========\n";

  return 0;
}
