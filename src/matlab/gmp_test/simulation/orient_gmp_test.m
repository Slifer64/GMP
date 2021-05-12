% function orient_gmp_test()

clc;
close all;
clear;

set_matlab_utils_path('../');

%% Load training data

load('orient_data.mat', 'Data');
Timed = Data.Time;
Qd_data = Data.Quat;
vRotd_data = Data.rotVel;
dvRotd_data = Data.rotAccel;

% fid = FileIO('orient_data.bin', bitor(FileIO.out, FileIO.trunc) );
% fid.write('Timed',Timed);
% fid.write('Qd_data',Qd_data);
% fid.write('vRotd_data',vRotd_data);
% fid.write('dvRotd_data',dvRotd_data);
% fid.close();

%% Write data to binary format
% fid = fopen('gmp_orient_train_data.bin','w');
% io_.write_mat(Timed, fid, true);
% io_.write_mat(Qd_data, fid, true);
% io_.write_mat(vRotd_data, fid, true);
% io_.write_mat(dvRotd_data, fid, true);
% fclose(fid);

Ts = Timed(2)-Timed(1);

simulateGMPo = @simulateGMPo_in_Cart_space; % simulateGMPo_in_'log/quat/Cart'_space

%% initialize and train GMP
train_method = 'LS';
N_kernels = 30;
kernels_std_scaling = 1;
gmp_o = GMPo(N_kernels, kernels_std_scaling);
tic
offline_train_mse = gmp_o.train(train_method, Timed/Timed(end), Qd_data);
offline_train_mse
toc

% traj_sc = TrajScale_Prop(n_dof);
% traj_sc = TrajScale_Rot_min();
traj_sc = TrajScale_Rot_wb();
traj_sc.setWorkBenchNormal([0; 0; 1]);
gmp_o.setScaleMethod(traj_sc);

% gmp_o.exportToFile('data/gmp_o_model.bin');
% gmp_o = GMPo.importFromFile('data/gmp_o_model.bin');

%% DMP simulation
disp('GMP simulation...');
tic
Qd0 = Qd_data(:,1);
Q0 = Qd0; %math_.quatProd(rotm2quat(rotz(57))',Qd0);
Qgd = Qd_data(:,end);
ks = diag([1.2 1.3 1.1]);
kt = 2;
e0 = ks*gmp_.quatLog(gmp_.quatDiff(Qgd,Qd0));
Qg = gmp_.quatProd(gmp_.quatExp(e0), Q0);
T = Timed(end) / kt;
dt = Ts;

% gmp_o.setScaleMethod(TrajScale.PROP_SCALE);

[Time, Q_data, vRot_data, dvRot_data] = simulateGMPo(gmp_o, Q0, Qg, T, dt);
toc

% Data.Time = Time;
% Data.Quat = Q_data;
% Data.rotVel = vRot_data;
% Data.rotAccel = dvRot_data;
% save('data/orient_train_data.mat', 'Data');

%% Plot results

Timed = Timed/kt;

T_sc = gmp_o.getScaling();

Pqd_data0 = zeros(3, size(Qd_data,2));
Pqd_data = zeros(3, size(Qd_data,2));
for j=1:size(Pqd_data,2)
    Pqd_data0(:,j) = GMPo.quat2q(Qd_data(:,j), Qd0);
    Pqd_data(:,j) = T_sc*Pqd_data0(:,j);
    Qd_data(:,j) = GMPo.q2quat(Pqd_data(:,j), Q0);
end

for j=1:size(vRotd_data,2)-1
    vRotd_data(:,j) = math_.quatLog(math_.quatDiff(Qd_data(:,j+1),Qd_data(:,j)))/Ts;
end
vRotd_data(:,j) = zeros(3,1);
for i=1:3, dvRotd_data(i,:)=[diff(vRotd_data(i,:)) 0]/Ts; end


Pq_data = zeros(3, size(Q_data,2));
for j=1:size(Pq_data,2)
    Pq_data(:,j) = GMPo.quat2q(Q_data(:,j), Q0);
end

line_width = 2.5;

figure('Position', [200 200 600 500]);
y_labels = {'$e_{q,x}$','$e_{q,y}$', '$e_{q,z}$'};
for i=1:3
   subplot(3,1,i);
   hold on;
   plot(Time, Pq_data(i,:), 'LineWidth', line_width, 'Color','blue');
   plot(Timed, Pqd_data(i,:), 'LineWidth', line_width, 'LineStyle',':', 'Color',[0.85 0.33 0.1 0.85]);
   plot(Timed, Pqd_data0(i,:), 'LineWidth', line_width, 'LineStyle','--', 'Color',[0.93 0.69 0.13 0.7]);
   ylabel(y_labels{i}, 'interpreter','latex', 'fontsize',20);
   axis tight;
   if (i==1), legend({'$gmp$', '$k_s * demo$', '$demo$'}, 'interpreter','latex', 'fontsize',16, 'Position',[0.7 0.78 0.27 0.15]); end
   if (i==1), title('Quaternion error: $e_q = log(Q * Q_0^{-1})$', 'interpreter','latex', 'fontsize',18); end
   if (i==3), xlabel('time [$s$]', 'interpreter','latex', 'fontsize',17); end
   hold off;
end

figure;
hold on;
plot3(Pq_data(1,:), Pq_data(2,:), Pq_data(3,:), 'LineWidth', line_width, 'LineStyle','-', 'Color','blue');
plot3(Pqd_data(1,:), Pqd_data(2,:), Pqd_data(3,:), 'LineWidth', line_width, 'LineStyle',':', 'Color',[0.85 0.33 0.1]);
plot3(Pqd_data0(1,:), Pqd_data0(2,:), Pqd_data0(3,:), 'LineWidth', line_width, 'LineStyle','--', 'Color',[0.93 0.69 0.13]);
legend({'$gmp$', '$k_s * demo$', '$demo$'}, 'interpreter','latex', 'fontsize',17);
hold off;


figure;
Q_labels = {'$\eta$','$\epsilon_1$', '$\epsilon_2$', '$\epsilon_3$'};
Qd_labels = {'$\eta_d$','$\epsilon_{d,1}$', '$\epsilon_{d,2}$', '$\epsilon_{d,3}$'};
for i=1:4
   subplot(4,1,i);
   hold on;
   plot(Time, Q_data(i,:), 'LineWidth', line_width);
   plot(Timed, Qd_data(i,:), 'LineWidth', line_width, 'LineStyle',':');
   legend({Q_labels{i}, Qd_labels{i}}, 'interpreter','latex', 'fontsize',15);
   if (i==1), title('Unit Quaternion', 'interpreter','latex', 'fontsize',17); end
   if (i==4), xlabel('time [$s$]', 'interpreter','latex', 'fontsize',15); end
   hold off;
end

return

figure;
vRot_labels = {'$\omega_x$','$\omega_y$', '$\omega_z$'};
vRotd_labels = {'$\omega_{d,x}$','$\omega_{d,y}$', '$\omega_{d,z}$'};
for i=1:3
   subplot(3,1,i);
   hold on;
   plot(Time, vRot_data(i,:), 'LineWidth', line_width);
   plot(Timed, vRotd_data(i,:), 'LineWidth', line_width, 'LineStyle',':');
   legend({vRot_labels{i}, vRotd_labels{i}}, 'interpreter','latex', 'fontsize',15);
   if (i==1), title('Rotational Velocity', 'interpreter','latex', 'fontsize',17); end
   if (i==3), xlabel('time [$s$]', 'interpreter','latex', 'fontsize',15); end
   hold off;
end

figure;
dvRot_labels = {'$\dot{\omega}_x$','$\dot{\omega}_y$', '$\dot{\omega}_z$'};
dvRotd_labels = {'$\dot{\omega}_{d,x}$','$\dot{\omega}_{d,y}$', '$\dot{\omega}_{d,z}$'};
for i=1:3
   subplot(3,1,i);
   hold on;
   plot(Time, dvRot_data(i,:), 'LineWidth', line_width);
   plot(Timed, dvRotd_data(i,:), 'LineWidth', line_width, 'LineStyle',':');
   legend({dvRot_labels{i}, dvRotd_labels{i}}, 'interpreter','latex', 'fontsize',15);
   if (i==1), title('Rotational Acceleration', 'interpreter','latex', 'fontsize',17); end
   if (i==3), xlabel('time [$s$]', 'interpreter','latex', 'fontsize',15); end
   hold off;
end

% Data.Time = Time;
% Data.Quat = Q_data;
% Data.rotVel = vRot_data;
% Data.rotAccel = dvRot_data;
% save('data/orient_train_data_overscale.mat', 'Data');


% end
