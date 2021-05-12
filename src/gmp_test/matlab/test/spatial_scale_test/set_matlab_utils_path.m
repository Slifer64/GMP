function set_matlab_utils_path()

path = strrep(mfilename('fullpath'), 'set_matlab_utils_path','');

addpath([path '../utils/']);
addpath([path '../../lib/dmp_lib/']);
addpath([path '../../lib/io_lib/']);
addpath([path '../../lib/math_lib/']);
addpath([path '../../lib/plot_lib/']);

import_dmp_lib();
import_io_lib();
import_math_lib();
import_plot_lib();

end
