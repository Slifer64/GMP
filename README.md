# GMP - Generalizable Movement Primitives Library.

GMP is library that provides an interface for DMP (Dynamic Movement Primitives). It is based on the DMP formulation presented in [this paper](https://arxiv.org/abs/2010.07708)

## Main Features
+ Encoding and generalization of a kinematic behaviour.
+ Online adaptation to via-points.
+ Off-line optimization under kinematic constraints.
+ Onl-ine MPC-like optimization under kinematic constraints.

## Supported Interfaces
+ matlab
+ c++

Both interfaces are implemented in a similar way employing the same classes and functions (except for some small exceptions due to the specific structure of each language).

## Dependencies
+ [armadillo](http://arma.sourceforge.net/)
It's best to have the latest version of armadillo installed. Older versions might cause strange runtime errors (that are not related to armadillo in any obvious way).
+ [OSQP library](https://github.com/osqp/osqp)
Binaries are already included in the *osqp_lib* folder which is ros package wrapper for the libary.
To update the library to the latest version, run:

```sh
  ./install/osqp_install.sh
```

## Examples
