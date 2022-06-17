# GMP - Generalizable Movement Primitives Library.

GMP is library that provides an interface for DMP (Dynamic Movement Primitives). It is based on the DMP formulation presented in [this paper](https://ieeexplore.ieee.org/abstract/document/9562059)

---

## Main Features
+ Encoding and generalization of a kinematic behaviour.
+ Online adaptation to via-points.
+ Off-line optimization under kinematic constraints.
+ Onl-ine MPC-like optimization under kinematic constraints.

---

## Supported Interfaces
+ matlab
+ c++

Both interfaces are implemented in a similar way employing the same classes and functions (except for some small exceptions due to the specific structure of each language).
The matlab libraries are located in `src/lib/matlab` and the c++ libraries in `src/lib/c++`.

---

## Dependencies for C++ interface
+ catkin
+ [armadillo](http://arma.sourceforge.net/)
  It's best to have the latest version of armadillo installed. Older versions might cause strange runtime errors (that are not related to armadillo in any obvious way).
+ [OSQP library](https://github.com/osqp/osqp)
  Binaries are already included in the *osqp_lib* folder which is ros package wrapper for the libary.
  To update the library to the latest version, run:

    ```sh
      ./install/osqp_install.sh
    ```

---

## Examples

Examples are provided in the package `src/gmp_test` both in matlab and c++.

### Setup
To run the matlab code, first setup your matlab search path. To do so, open a matlab,
browse to `src\gmp_test\src` and from the matlab command prompt type:
```matlab
  add_libraries_to_matlab_path
```
To run the c++ code, you have to build it first. To do so, open a terminal,
browse to the repository head and type:
 ```sh
  ./build.sh
```

### Run

+ #### Encoding and reproduction of a Cartesian position trajectroy:

    **matlab**: from matlab, browse to `src\gmp_test\src\simulation` and run:
    ```matlab
      nDoF_gmp_test
    ```
    **c++**: open a terminal, browse to the repository head and run:
    ```sh
      source devel/setup.bash
      roslaunch gmp_test nDoF_GMP_test.launch
    ```
    To plot the results, from matlab browse to `src\gmp_test\src\simulation` and run:
    ```matlab
      nDoF_gmp_plot
    ```
    You can also alter the simulation parameters by changing the ros parameters from `src\gmp_test\launch\nDoF_GMP_test.launch`
    
+ #### Encoding and reproduction of a Cartesian orientation trajectroy:
     **matlab**: browse from matlab to `src\gmp_test\src\simulation` and run:
    ```matlab
      orient_gmp_test
    ```
    **c++**: open a terminal, browse to the repository head and run:
    ```sh
      source devel/setup.bash
      roslaunch gmp_test orient_GMP_test.launch
    ```
    To plot the results, browse from matlab to `src\gmp_test\src\simulation` and run:
    ```matlab
      orient_gmp_plot
    ```
    You can also alter the simulation parameters by changing the ros parameters from `src\gmp_test\launch\orient_GMP_test.launch`
    
+ #### Adding via-points:

+ #### Enforcing constraints offline:

+ #### Enforcing constraints online:

---

### Miscellaneous

+ #### IO functionalities

    The package has additionally IO functionalities for reading/writing data from/to binary files, provided by the class `FileIO`. Each variable has a name-identifier assigned to it, which is used to write or read the variable. In general, it works like this for writing:
    ```matlab
    % data to write
    A = rand(4,5);
    v = rand(1,3);
    str = 'Hello World!';
    gmp = GMP(6, 30); % DMP with 6 DoFs and 30 kernels per DoF
    gmp_o = GMPo(40); % orientation DMP with 40 kernels
    
    % write to file
    fid = FileIO('my_data.bin', FileIO.out);
    fid.write('my_vector', v);
    fid.write('my_message', str);
    fid.write('some_matrix', A);
    gmp_.write(gmp, fid, 'my_gmp'); % use dedicated function for writing GMP
    gmp_.write(gmp_o, fid, 'o_gmp');
    ```
    and to read:
    ```matlab
    % read from file (in any order you like)
    fid = FileIO('my_data.bin', FileIO.in);
    
    gmp = GMP(); % create empty object
    gmp_.read(gmp, fid, 'my_gmp'); % use dedicated function for reading GMP

    A = fid.read('some_matrix');

    gmp_o = GMPo(); % create empty object
    gmp_.read(gmp_o, fid, 'o_gmp');
    
    v = fid.read('my_vector');
    str = fid.read('my_message');
    ```
    To see more details on how it works, check the examples `io_example.m` and `gmp_io_example.m` in the folder `src\gmp_test\src\simulation`. The API is almost identical for c++ and compatible with the matlab API. Therefore you can read from c++, data written in matlab and vice versa.  



---

## Citing GMP

If you are using GMP for your work, you are encouraged to

- Cite the related papers:
  ```
  @INPROCEEDINGS{9562059,  
  author={Sidiropoulos, Antonis and Doulgeri, Zoe},  
  booktitle={2021 IEEE International Conference on Robotics and Automation (ICRA)},   
  title={A Reversible Dynamic Movement Primitive formulation},   
  year={2021},  
  volume={},  
  number={},  
  pages={3147-3153},  
  doi={10.1109/ICRA48506.2021.9562059}}
  ```
- Put a star on this repository.
