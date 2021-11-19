
matlab_lib_path = [pwd '/../../lib/matlab/'];

%% add to seach path for the current matlab session
addpath(matlab_lib_path);

%% add it permantly
% savepath(matlab_lib_path);

%% restore path, if you add it permantly and want to remove it
% rmpath(matlab_lib_path);
% savepath;