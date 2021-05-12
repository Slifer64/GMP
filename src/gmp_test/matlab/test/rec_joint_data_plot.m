clc;
close all;
clear;

set_matlab_utils_path();

% fid = FileIO('joints_rec_data.bin', FileIO.in);
% fid.printHeader();
% jpos_data = fid.read('joints');
% plot(jpos_data')
% return

fid = FileIO('data/joints_rec_data.bin', FileIO.in);
fid.printHeader();

Time = fid.read('time');
joints_data = fid.read('joints');

time_for = [0.95; ];

figure;
plot(Time, joints_data(1,:));

return


fid = FileIO('data/train_data.bin', FileIO.in);
% fid = FileIO('data/sim_data.bin', FileIO.in);

fid.printHeader();

joint_pos_data = fid.read('joint_pos_data');
Pd_data = fid.read('Pd_data');
Qd_data = fid.read('Qd_data');


%% -----------------------------------------------

j0 = joint_pos_data(:,1);
fprintf('j0: [');
for i=1:(length(j0)-1)
    fprintf('%.3f, ', j0(i));
end
fprintf('%.3f]\n',j0(end));

%% -----------------------------------------------
target_pose = [Pd_data(:,end); Qd_data(:,end)];
fprintf('target_pose: [');
for i=1:(length(target_pose)-1)
    fprintf('%.3f; ', target_pose(i));
end
fprintf('%.3f]\n',target_pose(end));
