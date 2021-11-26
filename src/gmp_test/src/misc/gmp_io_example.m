%% Example for reading/writing DMP from/to file
% (that can also contain other data as well)

clc;
close all;
clear;

import_io_lib();
import_gmp_lib();

filename = 'my_data.bin';

write_data_example(filename);

read_data_example(filename);

% cleanup
delete(filename);


%% =========== Write example =============
function write_data_example(filename)
    
    disp('======== Write example =========');
    
    %% data to write
    M = [1.1 2.2 3.3; 4.4 5.5 6.6; -7.7 8.8 -9.9];
    v = [1.1; 2.2; 3.3; 4.4; 5.5];
    str = 'GMP io example';
    
    gmp = GMP(3, 30); % create GMP with 3 DoFs and 30 kernels (per DoF)
    gmp_o = GMPo(25); % create orientation GMP with 25 kernels
    
    % Train the GMPs, do other stuff...

    %% Create a @FileIO object with open-mode 'out|trunc' for writing to it.
    % If a previous file with the same name existed, it will be deleted due
    % to 'trunc'.
    fid = FileIO(filename, bitor(FileIO.out,FileIO.trunc));
    
    % If we want to keep the contents of a previous file with same name then: 
    % fid = FileIO('my_data.bin', FileIO.out);
  
    %% write the variables in any order you like (doesn't matter)
    % by assigning a name-identifier to each one
    
    fprintf('Writing data...');
    
    fid.write('some_matrix', M);
    
    % The GMP have special structure, so use these dedicated functions for
    % writing them to FileIO object. The use of a name-identifier is in
    % general optional. However, it is highly recommended in cases where 
    % you write the gmp to a file which contains other variables too, in 
    % order to avoid name conficts.
    gmp_.write(gmp, fid, 'my_gmp'); % add name-identifier to avoid possible name conficts with other variables
    gmp_.write(gmp_o, fid, 'my_orient_gmp'); % same here...
    
    fid.write('a_vector', v);
    
    fprintf('[SUCCESS]\n');

end

%% =========== Read example =============
function read_data_example(filename)

    disp('======== Read example =========');
    
    %% Create a @FileIO object with open-mode 'in' for reading from it.
    fid = FileIO(filename, FileIO.in);
    
    disp('File contents info:');
    
    % To view the name and types of the variables
    % fid.printHeader();
    
    fprintf('Reading data...');

    % read the variables in any order you like
    % providing the name-identifier of each one
    
    v = fid.read('a_vector');
    M = fid.read('some_matrix');
    
    %% read the GMPs
    gmp_o = GMPo(); % create an empty object first
    gmp_.read(gmp_o, fid,'my_orient_gmp'); % and then read it
    
    gmp = GMP();
    gmp_.read(gmp, fid,'my_gmp');
    

    fprintf('[SUCCESS]\n');

end
