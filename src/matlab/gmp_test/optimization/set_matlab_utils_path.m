function set_matlab_utils_path(prefix)

if (nargin < 1), prefix=''; end

path = strrep(mfilename('fullpath'), 'set_matlab_utils_path','');

addpath([path '/' prefix 'utils/']);
addpath([path '/' prefix '../lib/gmp_lib/']);
addpath([path '/' prefix '../lib/io_lib/']);
addpath([path '/' prefix '../lib/math_lib/']);
addpath([path '/' prefix '../lib/plot_lib/']);

import_gmp_lib();
import_io_lib();
import_math_lib();
import_plot_lib();

end
