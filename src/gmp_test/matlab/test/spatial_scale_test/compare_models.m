clc;
close all;
clear;

format compact;

set_matlab_utils_path;

%% =====  load train data  =====
fid = FileIO('../data/ur_train_data.bin', FileIO.in);
% fid.printHeader();
Time = fid.read('Timed');
Pd_data = fid.read('Pd_data');
dPd_data = fid.read('dPd_data');
Qd_data = fid.read('Qd_data');
% vRotd_data = fid.read('vRotd_data');

%% =====  train model  =====
target_model = TargetModel();
target_model_train_err = target_model.train(Time, Pd_data, Qd_data);

std_model = StdModel();
std_model_train_err = std_model.train(Time, Pd_data, Qd_data);

rot_model = RotModel();
rot_model_train_err = rot_model.train(Time, Pd_data, Qd_data);

%% =====  simulate  =====
P0 = Pd_data(:,1);
Pg = Pd_data(:,end);
Q0 = Qd_data(:,1);
Qg = Qd_data(:,end);
Pg = 1.5*(Pg-P0) + P0;
Qg = math_.quatProd(rotm2quat(rotz(90))', Qg);
tau = Time(end) - Time(1);
x = Time/tau;

target_model.simulate(P0, Q0, Pg, Qg, tau, x);
std_model.simulate(P0, Q0, Pg, Qg, tau, x);
rot_model.simulate(P0, Q0, Pg, Qg, tau, x);

% ax_legend = {{'$x$','$x_d$'}, {'$y$','$y_d$'}, {'$z$','$z_d$'}};
% fig = figure;
% fig.Position(3:4) = [1444 865];
% for i=1:3
%     subplot(2,3,i); hold on;
%     plot(Time, P_data(i,:), 'LineWidth',2, 'Color','blue');
%     plot(Time, Pd_data(i,:), 'LineWidth',2, 'Color','magenta', 'LineStyle','--');
%     legend(ax_legend{i}, 'interpreter','latex', 'fontsize',15);
%     if (i==1), ylabel('World pos [$m$]', 'interpreter','latex', 'fontsize',15); end
%         
%     subplot(2,3,i+3); hold on;
%     plot(Time, P_t_data(i,:), 'LineWidth',2, 'Color','blue');
%     plot(Time, Pd_t_data(i,:), 'LineWidth',2, 'Color','magenta', 'LineStyle','--');
%     legend(ax_legend{i}, 'interpreter','latex', 'fontsize',15);
%     xlabel('Time [$s$]', 'interpreter','latex', 'fontsize',15);
%     if (i==1), ylabel('Target pos [$m$]', 'interpreter','latex', 'fontsize',15); end    
% end

%% ===================================================================
% Position

figure; 
ax = axes(); hold on;
plot3(target_model.P_data(1,:), target_model.P_data(2,:), target_model.P_data(3,:), 'LineWidth',2, 'LineStyle','-', 'Color',[0.85 0.33 0.1]);
plot3(rot_model.P_data(1,:), rot_model.P_data(2,:), rot_model.P_data(3,:), 'LineWidth',2, 'LineStyle','--', 'Color','magenta');
plot3(std_model.P_data(1,:), std_model.P_data(2,:), std_model.P_data(3,:), 'LineWidth',2, 'LineStyle','-.', 'Color','blue');
plot3(Pd_data(1,:), Pd_data(2,:), Pd_data(3,:), 'LineWidth',2, 'LineStyle','--', 'Color','cyan');
scatter3(P0(1), P0(2), P0(3), 'LineWidth',4, 'SizeData',100, 'Marker','o' , 'MarkerEdgeColor','green');
scatter3(Pg(1), Pg(2), Pg(3), 'LineWidth',4, 'SizeData',100, 'Marker','x' , 'MarkerEdgeColor','red');
plotFrame(ax, Pg, Qg, 0.1);
plotFrame(ax, P0, Q0, 0.1);
plotFrame(ax, Pd_data(:,end), Qd_data(:,end), 0.1);
% plotFrames(ax, target_model.P_data, target_model.Q_data, 0.1);
legend({'$p_t$','$p_{std}$', '$p_{rot}$', '$p_d$', '$p_0$', '$p_g$'}, 'interpreter','latex', 'fontsize',15);
plotAnimatedFrame(ax, rot_model.P_data, rot_model.Q_data, 10, 1e-3, 0.1);
plotAnimatedFrame(ax, target_model.P_data, target_model.Q_data, 10, 1e-3, 0.1);

return

%% ===================================================================
% Orientation

qg = GMPo.quat2q(Qg, Q0);
q0 = GMPo.quat2q(Q0, Q0);
n_data = size(Qd_data,2);
qd_data = zeros(3, n_data);
for j=1:n_data, qd_data(:,j) = GMPo.quat2q(Qd_data(:,j), Q0); end

figure; 
ax = axes(); hold on;
plot3(target_model.q_data(1,:), target_model.q_data(2,:), target_model.q_data(3,:), 'LineWidth',2, 'LineStyle','-', 'Color',[0.85 0.33 0.1]);
plot3(rot_model.q_data(1,:), rot_model.q_data(2,:), rot_model.q_data(3,:), 'LineWidth',2, 'LineStyle','--', 'Color','magenta');
plot3(std_model.q_data(1,:), std_model.q_data(2,:), std_model.q_data(3,:), 'LineWidth',2, 'LineStyle','-.', 'Color','blue');
plot3(qd_data(1,:), qd_data(2,:), qd_data(3,:), 'LineWidth',2, 'LineStyle','--', 'Color','cyan');
scatter3(q0(1), q0(2), q0(3), 'LineWidth',4, 'SizeData',100, 'Marker','o' , 'MarkerEdgeColor','green');
scatter3(qg(1), qg(2), qg(3), 'LineWidth',4, 'SizeData',100, 'Marker','x' , 'MarkerEdgeColor','red');
plotFrame(ax, qg, Qg, 0.1);
plotFrame(ax, q0, Q0, 0.1);
plotFrame(ax, qd_data(:,end), Qd_data(:,end), 0.1);
legend({'$p_t$','$p_{std}$', '$p_{rot}$', '$p_d$', '$p_0$', '$p_g$'}, 'interpreter','latex', 'fontsize',15);

%% ===================================================================


function plotFrames(ax, P_data, Q_data, frame_scale)

n_data = size(Q_data,2);
ind = floor([0.2 0.3 0.4 0.5 0.6 0.7] * n_data);

for j=1:length(ind)
    plotFrame(ax, P_data(:,ind(j)), Q_data(:,ind(j)), frame_scale);
end

end


