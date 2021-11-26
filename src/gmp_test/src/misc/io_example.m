%% Example on how @FileIO works.
% You can use it to read/write variables from/to a binary file.
% A name-identifier is assigned to each variable, which can be used to
% read/write the corresponding variable.
% The supported types are: scalars, matrices, character arrays (strings)
% Open modes: 'in', 'out', 'trunc'
% Duplicate name-identifiers are disallowed.

clc;
close all;
clear;

import_io_lib();

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
    col_vec = [1.1; 2.2; 3.3; 4.4; 5.5];
    row_vec = -[1.1, 2.2, 3.3, 4.4, 5.5];
    str = 'Hello world';
    d = 3.14159; % double
    i = int32(2); % int
    long_i = int64(10000000000); % long int
    Mui = uint32( [1 2; 3 4] ); % matrix of type int
    

    %% Create a @FileIO object with open-mode 'out|trunc' for writing to it.
    % If a previous file with the same name existed, it will be deleted due
    % to 'trunc'.
    fid = FileIO(filename, bitor(FileIO.out,FileIO.trunc));
    
    % If we want to keep the contents of a previous file with same name then: 
    % fid = FileIO('my_data.bin', FileIO.out);

    %% write the variables in any order you like (doesn't matter)
    % by assigning a name-identifier to each one
    
    fprintf('Writing data...');
    
    fid.write('a_vec', col_vec);
    fid.write('my_matrix', M);
    fid.write('message', str);
    fid.write('a_row_vec', row_vec);
    fid.write('some_double', d);
    fid.write('some_integer', i);
    fid.write('some_long_integer', long_i);
    fid.write('mat_of_uint', Mui);
    
    fprintf('[SUCCESS]\n');
    
    disp('File contents info:');
    
    %% you can also view the variables that have been written so far and
    % their type by printing the header:
    fid.printHeader();

    %% (optionally) close the file
    fid.close();
    % this operation is optional. It will close automatically when
    % the @FileIO object goes out of scope.

end

%% =========== Read example =============
function read_data_example(filename)

    disp('======== Read example =========');
    
    %% Create a @FileIO object with open-mode 'in' for reading from it.
    fid = FileIO(filename, FileIO.in);
    
    disp('File contents info:');
    
    %% To view the name and types of the variables contained in
    % 'my_data.bin':
    fid.printHeader();

    %% read the variables in any order you like
    % providing the name-identifier of each one
    
    fprintf('Reading data...');
    
    col_vec = fid.read('a_vec');
    M = fid.read('my_matrix');
    i = fid.read('some_integer');
    row_vec = fid.read('a_row_vec');
    d = fid.read('some_double');
    Mui = fid.read('mat_of_uint');
    long_i = fid.read('some_long_integer');
    str = fid.read('message');
    
    %% alternatively you can get all the above data as a struct:
    % data = fid.readAll();
    % data.a_vec
    % data.my_matrix
    % ...
    
    fprintf('[SUCCESS]\n');

    %% (optionally) close the file
    fid.close();
    % this operation is optional. It will close automatically when
    % the @FileIO object goes out of scope.
    
    %% Validation:
    err_tol = 1e-16;
    
    expected_M = [1.1 2.2 3.3; 4.4 5.5 6.6; -7.7 8.8 -9.9];
    if (norm(M-expected_M,'fro') > err_tol), error('M is different!'); end

    expected_col_vec = [1.1; 2.2; 3.3; 4.4; 5.5];
    if (norm(col_vec-expected_col_vec) > err_tol), error('col_vec is different!'); end
    
    expected_row_vec = -[1.1, 2.2, 3.3, 4.4, 5.5];
    if (norm(row_vec-expected_row_vec) > err_tol), error('row_vec is different!'); end
    
    expected_str = 'Hello world';
    if (norm(str-expected_str) > err_tol), error('str is different!'); end
    
    expected_d = 3.14159; % double
    if (norm(M-expected_M,'fro') > err_tol), error('M is different!'); end
    
    expected_i = int32(2); % int
    if (abs(d-expected_d) > err_tol), error('d is different!'); end
    
    expected_long_i = int64(10000000000); % long int
    if (abs(long_i-expected_long_i) > err_tol), error('long_i is different!'); end
    
    expected_Mui = uint32( [1 2; 3 4] ); % matrix of type int
    if (norm(double(Mui-expected_Mui),'fro') > err_tol), error('Mui is different!'); end
    
    disp('Validation success!');
    
end
