clc;
close all;
clear;

set_matlab_utils_path('../');

%% Load training data

load('pos_data.mat', 'Data');
Timed = Data.Time;
P_data = Data.Pos;

Ts = Timed(2)-Timed(1);


%% initialize and train GMP
train_method = 'LS';
N_kernels = 30;
kernels_std_scaling = 1;
gmp = GMP(3, N_kernels, kernels_std_scaling);
tic
x = Timed/Timed(end);
offline_train_mse = gmp.train(train_method, x, P_data);
offline_train_mse
toc

% traj_sc = TrajScale_Prop(3);
% traj_sc = TrajScale_Rot_min();
% traj_sc = TrajScale_Rot_wb();
% traj_sc.setWorkBenchNormal([0; 0; 1]);
% gmp_o.setScaleMethod(traj_sc);

n_data = length(x);

%% Construct new orientation trajectory
P0 = P_data(:,1);
Pg = P_data(:,end);

gmp.setY0(P0);
gmp.setGoal(Pg);

% gmp.setY0(P0);

Pd_new = zeros(3, n_data);
dP_new = zeros(3, n_data);
ddP_new = zeros(3, n_data);

for j=1:n_data
    Pd_new(:,j) = P_data(:,j);% + 0.3*exp(-((0.4-x(j))/0.1)^2) - 0.2*exp(-((0.7-x(j))/0.06)^2);
end

for i=1:3
    dP_new(i,:) = [diff(Pd_new(i,:)) 0]/Ts;
    ddP_new(i,:)=[diff(dP_new(i,:)) 0]/Ts;
end

%% Update GMP

x_dot = 1/Timed(end);
x_ddot = 0;

gmp_new = gmp.deepCopy();
gmp_up = GMP_Update(gmp_new);
gmp_up.enableSigmawUpdate(false);
gmp_up.setMsrNoiseVar(1e-4);

up_ind = [2324];
for j=1:length(up_ind)
    gmp_up.updatePos(x(up_ind(j)), Pd_new(:,up_ind(j)));
end


% for j=1:10:n_data
%     j/n_data
%     %gmp_up.updateQuat(x, Qd_new(:,j));
%     gmp_up.updatePos(x, qd_new_data(:,j));
% %    gmp_up.updateRotVel(x, x_dot, vRot_new(:,j), Q_new(:,j) );
% %    gmp_up.updateRotAccel(x, x_dot, dvRot_new(i,:), vRot_new(:,j), Q_new(:,j) );
% end

P_data = zeros(3, n_data);
P_new_data = zeros(3, n_data);
for j=1:n_data
    P_data(:,j) = gmp.getYd(x(j));
    P_new_data(:,j) = gmp_new.getYd(x(j));
end


figure;
for i=1:3
   subplot(3,1,i); hold on;
   plot(Timed, Pd_new(i,:), 'LineWidth',2, 'Color','blue', 'LineStyle','-');
   plot(Timed, P_data(i,:), 'LineWidth',2, 'Color',[0.85 0.33 0.1], 'LineStyle','--');
   plot(Timed, P_new_data(i,:), 'LineWidth',2, 'Color','magenta', 'LineStyle',':');
   plot(Timed(up_ind), Pd_new(i,up_ind), 'LineStyle','none', 'Marker','*', 'MarkerSize',10, 'Color','red', 'LineWidth',2, 'HandleVisibility','off');
   if (i==1), legend({'$P_{d,new}$', '$P$', '$P_{new}$'}, 'interpreter','latex', 'fontsize',15); end
   if (i==4), xlabel('time [$s$]', 'interpreter','latex', 'fontsize',15); end
end


