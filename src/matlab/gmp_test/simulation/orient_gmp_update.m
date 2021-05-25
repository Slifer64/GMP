% function orient_gmp_test()

clc;
close all;
clear;

set_matlab_utils_path('../');

%% Load training data

load('orient_data.mat', 'Data');
Timed = Data.Time;
Qd_data = Data.Quat;
% vRot_new = Data.rotVel;
% dvRot_new = Data.rotAccel;

Ts = Timed(2)-Timed(1);


%% initialize and train GMP
train_method = 'LS';
N_kernels = 30;
kernels_std_scaling = 1;
gmp_o = GMPo(N_kernels, kernels_std_scaling);
tic
x = Timed/Timed(end);
offline_train_mse = gmp_o.train(train_method, x, Qd_data);
offline_train_mse
toc

% traj_sc = TrajScale_Prop(3);
% traj_sc = TrajScale_Rot_min();
% traj_sc = TrajScale_Rot_wb();
% traj_sc.setWorkBenchNormal([0; 0; 1]);
% gmp_o.setScaleMethod(traj_sc);

n_data = length(x);

%% Construct new orientation trajectory
Q0 = Qd_data(:,1);

% gmp_o.setQ0(Q0);

qd_data = zeros(3, n_data);

Qd_new = zeros(4, n_data);
qd_new_data = zeros(3, n_data);
vRotd_new = zeros(3, n_data);
dvRotd_new = zeros(3, n_data);

for j=1:n_data
    qd_data(:,j) = GMPo.quat2q(Qd_data(:,j), Q0);
    qd_new_data(:,j) = qd_data(:,j) + 0.3*exp(-((0.4-x(j))/0.1)^2) - 0.2*exp(-((0.7-x(j))/0.06)^2);
    Qd_new(:,j) = GMPo.q2quat(qd_new_data(:,j), Q0);
end

for j=1:n_data-1
    vRotd_new(:,j) = math_.quatLog( math_.quatDiff(Qd_data(:,j+1),Qd_data(:,j)) )/Ts;
end
for i=1:3, dvRotd_new(i,:)=[diff(vRotd_new(i,:)) 0]/Ts; end


%% Update GMP

x_dot = 1/Timed(end);
x_ddot = 0;

gmp_o_new = gmp_o.deepCopy();
gmp_up = GMPo_Update(gmp_o_new);
gmp_up.enableSigmawUpdate(true);
gmp_up.setMsrNoiseVar(1e-4);

rng(0);
up_ind = [];
for j=1:n_data
    if ( norm(qd_new_data(:,j)-qd_data(:,j) ) > 0.05  && (rand() > 0.99))
        up_ind = [up_ind j];
    end
end
x_up = x(up_ind);

% for j=1:length(up_ind)
%     gmp_up.updatePos(x_up(j), qd_new_data(:,up_ind(j)));
% end


for j=1:length(up_ind)
    k = up_ind(j);
    gmp_up.updateQuat(x_up(j), Qd_new(:,k));
%     gmp_up.updatePos(x_up(j), qd_new_data(:,k));
%     gmp_up.updateRotVel(x_up(j), x_dot, vRotd_new(:,k), Qd_new(:,k) );
%     gmp_up.updateRotAccel(x_up(j), x_dot, x_ddot, dvRotd_new(:,k), vRotd_new(:,k), Qd_new(:,k) );
end

q_data = zeros(3, n_data);
Q_data = zeros(4, n_data);

q_new_data = zeros(3, n_data);
Q_new_data = zeros(4, n_data);
for j=1:n_data
    Q_data(:,j) = gmp_o.getQd(x(j));
    q_data(:,j) = gmp_o.getYd(x(j));
    
    Q_new_data(:,j) = gmp_o_new.getQd(x(j));
    q_new_data(:,j) = gmp_o_new.getYd(x(j));
end


figure;
for i=1:3
   subplot(3,1,i); hold on;
   plot(Timed, qd_new_data(i,:), 'LineWidth',2, 'Color','blue', 'LineStyle','-');
   plot(Timed, q_data(i,:), 'LineWidth',2, 'Color',[0.85 0.33 0.1], 'LineStyle','--');
   plot(Timed, q_new_data(i,:), 'LineWidth',2, 'Color','magenta', 'LineStyle',':');
   plot(Timed(up_ind), qd_new_data(i,up_ind), 'LineStyle','none', 'Marker','*', 'MarkerSize',10, 'Color','red', 'LineWidth',2, 'HandleVisibility','off');
   if (i==1), legend({'$q_{d,new}$', '$q$', '$q_{new}$'}, 'interpreter','latex', 'fontsize',15); end
   if (i==4), xlabel('time [$s$]', 'interpreter','latex', 'fontsize',15); end
end

% figure;
% for i=1:4
%    subplot(4,1,i); hold on;
%    plot(Timed, Qd_new(i,:), 'LineWidth',2, 'Color','blue', 'LineStyle','-');
%    plot(Timed, Q_data(i,:), 'LineWidth',2, 'Color',[0.85 0.33 0.1], 'LineStyle','--');
%    plot(Timed, Q_new_data(i,:), 'LineWidth',2, 'Color','magenta', 'LineStyle',':');
%    if (i==1), legend({'$Q_{d,new}$', '$Q$', '$Q_{new}$'}, 'interpreter','latex', 'fontsize',15); end
%    if (i==4), xlabel('time [$s$]', 'interpreter','latex', 'fontsize',15); end
% end

