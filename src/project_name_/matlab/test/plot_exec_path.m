clc;
close all;
clear;

set_matlab_utils_path();

fid = FileIO('data/exec_data.bin', FileIO.in);

fid.printHeader();

Time = fid.read('Time');
x_data = fid.read('x_data');
P_data = fid.read('P_data');
dP_data = fid.read('dP_data');
Q_data = fid.read('Q_data');
vRot_data = fid.read('vRot_data');
Fext_data = fid.read('Fext_data');
mf_ind = fid.read('motion_finish_ind');
mf_ind = mf_ind + 1;
mf_ind = [1 mf_ind length(Time)];

n_targets = ( length(mf_ind) - 1 ) / 2;

Q0 = Q_data(:,1);

Mt_data = cell(n_targets, 1);
k = 1;
for i=1:n_targets
    Mt_data{i} = struct('for',[], 'rev',[]);
    
    Mt_data{i}.for = struct('Time',[], 'x_data',[], 'Pos',[], 'Quat',[], 'qlog',[], 'Fext',[]);
    i1 = mf_ind(k);
    i2 = mf_ind(k+1)-1;
    Mt_data{i}.for.Time = Time(i1:i2);
    Mt_data{i}.for.x_data = x_data(i1:i2);
    Mt_data{i}.for.Pos = P_data(:, i1:i2);
    Mt_data{i}.for.Quat = Q_data(:, i1:i2);
    Quat = Mt_data{i}.for.Quat;
    qlog = zeros(3, size(Quat,2));
    for j=1:size(qlog,2), qlog(:,j) = math_.quatLog(math_.quatDiff(Quat(:,j),Q0)); end
    Mt_data{i}.for.qlog = qlog;
    Mt_data{i}.for.Fext = Fext_data(:, i1:i2);
    k = k + 1;
    
    Mt_data{i}.rev = struct('Time',[], 'x_data',[], 'Pos',[], 'Quat',[], 'qlog',[], 'Fext',[]);
    i1 = mf_ind(k);
    i2 = mf_ind(k+1)-1;
    Mt_data{i}.rev.Time = Time(mf_ind(k):mf_ind(k+1)-1);
    Mt_data{i}.rev.x_data = x_data(i1:i2);
    Mt_data{i}.rev.Pos = P_data(:, mf_ind(k):mf_ind(k+1)-1);
    Mt_data{i}.rev.Quat = Q_data(:, mf_ind(k):mf_ind(k+1)-1);
    Quat = Mt_data{i}.rev.Quat;
    qlog = zeros(3, size(Quat,2));
    for j=1:size(qlog,2), qlog(:,j) = math_.quatLog(math_.quatDiff(Quat(:,j),Q0)); end
    Mt_data{i}.rev.qlog = qlog;
    Mt_data{i}.rev.Fext = Fext_data(:, i1:i2);
    k = k + 1;
end


figure;
ax = plot_3Dpath_with_orientFrames(Mt_data{1}.for.Pos, Mt_data{1}.for.Quat, 'numberOfFrames',6, 'LineWidth',5, 'frameLineWidth',3, 'frameScale',0.3);
ax.FontSize = 12;
% grid on;


%% ==============================================================

