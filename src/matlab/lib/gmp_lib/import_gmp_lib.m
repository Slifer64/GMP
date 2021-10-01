
gmp_lib_path = strrep(mfilename('fullpath'), 'import_gmp_lib','');

addpath(gmp_lib_path);
addpath([gmp_lib_path '/GMP/']);
addpath([gmp_lib_path '/GMPo/']);
addpath([gmp_lib_path '/TrajScale/']);

addpath([gmp_lib_path '/deps/QP_Goldfarb_Idnani/']);
addpath([gmp_lib_path '/deps/OSQP/']);
