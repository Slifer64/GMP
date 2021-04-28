
math_lib_path = mfilename('fullpath');
math_lib_path = strrep(math_lib_path, 'import_math_lib','');

addpath(math_lib_path);
% addpath([math_lib_path '/quaternion/']);
