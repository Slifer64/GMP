#include <iostream>

#include <gmp_lib/GMP/GMP_MPC.h>

using namespace as64_;

int main(int argc, char **argv)
{
    gmp_::GMP::Ptr gmp(new gmp_::GMP());
    gmp_::GMP_MPC gmp_mpc( gmp.get(), 10, 0.02, 30, 1.5, {0,0,0}, {1e6,100,20});

    return 0;
}