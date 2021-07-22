clc;
close all;
clear;

addpath('../../../../matlab/lib/gmp_lib/');
import_gmp_lib();

addpath('../../../../matlab/lib/io_lib/');
import_io_lib();

%% Load training data
load('train_data_pos_constr.mat','Timed','Pd_data','dPd_data','ddPd_data');
% data = FileIO('train_data_pos_constr.bin', FileIO.in).readAll();
% Timed = data.Timed;
% Pd_data = data.Pd_data;

Ts = Timed(2)-Timed(1);

% dPd_data = zeros(size(Pd_data));
% ddPd_data = zeros(size(Pd_data));
% for i=1:size(Pd_data,1)
%     dPd_data(i,:) = [diff(Pd_data(i,:)) 0]/Ts;
%     ddPd_data(i,:) = [diff(dPd_data(i,:)) 0]/Ts;
% end


%% initialize and train GMP
train_method = 'LS';
N_kernels = 30;
kernels_std_scaling = 1.5;
n_dof = size(Pd_data,1);
gmp = GMP(n_dof, N_kernels, kernels_std_scaling);
tic
offline_train_mse = gmp.train(train_method, Timed/Timed(end), Pd_data);
offline_train_mse
toc

% gmp.setScaleMethod(TrajScale.ROT_MIN_SCALE);

taud = Timed(end);
yd0 = gmp.getYd(0); %Pd_data(1);
ygd = gmp.getYd(1); %Pd_data(end);

kt = 2.0; % temporal scaling
ks = diag([1 1 1]); % spatial scaling
tau = taud/kt;
y0 = yd0 + 0;
% yg = ks*(ygd - yd0) + y0;
yg = ygd + [0.1; -0.1; 0.2];

% ks = (yg - y0)/(ygd - yd0);
% ks

% gmp.setY0(y0);
% gmp.setGoal(yg);

% [Timed, Pd_data, dPd_data, ddPd_data] = getGMPTrajectory(gmp, tau, y0, yg);
% save('train_data_pos_constr.mat','Timed','Pd_data','dPd_data','ddPd_data');
% return;

gmp.setScaleMethod(TrajScale_Prop(3));
[Time, P_data, dP_data, ddP_data] = getGMPTrajectory(gmp, tau, y0, yg);

gmp.setScaleMethod(TrajScale_Rot_wb());
[Time, P_rot_wb_data, dP_rot_wb_data, ddP_rot_wb_data] = getGMPTrajectory(gmp, tau, y0, yg);

% gmp.setScaleMethod(TrajScale_Rot_min());
% [Time, P_rot_min_data, dP_rot_min_data, ddP_rot_min_data] = getGMPTrajectory(gmp, tau, y0, yg);

%% calculate scaled demo trajectory
T_sc = gmp.getScaling();
Time2 = Timed / kt;
P2_data = T_sc*(Pd_data - Pd_data(:,1)) + y0;
dP2_data = kt*T_sc*dPd_data;
ddP2_data = (kt^2)*T_sc*ddPd_data;

%% plot results

% plot 3D path
fig = figure;
fig.Position(3:4) = [815 716];
hold on;
plot3(y0(1), y0(2), y0(3), 'LineWidth', 4, 'LineStyle','none', 'Color','green','Marker','o', 'MarkerSize',8);
plot3(yg(1), yg(2), yg(3), 'LineWidth', 4, 'LineStyle','none', 'Color','red','Marker','x', 'MarkerSize',8);
plot3(ygd(1), ygd(2), ygd(3), 'LineWidth', 4, 'LineStyle','none', 'Color','magenta','Marker','x', 'MarkerSize',8);
plot3(P_data(1,:), P_data(2,:), P_data(3,:), 'LineWidth', 2, 'LineStyle','-', 'Color','blue');
plot3(P_rot_wb_data(1,:), P_rot_wb_data(2,:), P_rot_wb_data(3,:), 'LineWidth', 2, 'LineStyle','-', 'Color','cyan');
% plot3(P_rot_min_data(1,:), P_rot_min_data(2,:), P_rot_min_data(3,:), 'LineWidth', 2, 'LineStyle','-', 'Color','green');
plot3(Pd_data(1,:), Pd_data(2,:), Pd_data(3,:), 'LineWidth', 2, 'LineStyle',':', 'Color','magenta');
plot3([ygd(1) yg(1)], [ygd(2) yg(2)], [ygd(3) yg(3)], 'LineWidth', 1, 'LineStyle','--', 'Color','magenta');
legend({'$p_0$','$g$','$g_d$','$prop$','$rot-wb$','$demo$'}, 'interpreter','latex', 'fontsize',17);
view(170.5 ,5.05);
grid on;
hold off;

for i=1:n_dof
    figure;
    ax = cell(3,1);

    ax{1} = subplot(3,1,1);
    hold on;
    plot(Time, P_data(i,:), 'LineWidth',2.0, 'Color', 'blue');
    plot(Time, P_rot_wb_data(i,:), 'LineWidth',2.0, 'Color', 'cyan');
    plot(Time2, P2_data(i,:), 'LineWidth',2.0, 'Color', 'magenta', 'LineStyle',':');
    ylabel('pos [$m$]', 'interpreter','latex', 'fontsize',15);
    title(['temporal scale: $' num2str(kt) '$     ,     spatial scale: $' num2str(det(ks)) '$'], 'interpreter','latex', 'fontsize',18);
    scatter(nan, nan, 'MarkerEdgeColor',[0 0.7 0], 'Marker','o', 'MarkerEdgeAlpha',0.6, 'LineWidth',4, 'SizeData', 100); % dummy plot for legend
    scatter(nan, nan, 'MarkerEdgeColor',[1 0 0], 'Marker','^', 'MarkerEdgeAlpha',0.3, 'LineWidth',2, 'SizeData', 100); % dummy plot for legend
    scatter(nan, nan, 'MarkerEdgeColor',[1 0 0], 'Marker','v', 'MarkerEdgeAlpha',0.3, 'LineWidth',2, 'SizeData', 100); % dummy plot for legend
    % legend({'GMP','ref','eq-constr','low-bound','upper-bound'}, 'interpreter','latex', 'fontsize',15);
    axis tight;
    hold off;

    ax{2} = subplot(3,1,2);
    hold on;
    plot(Time, dP_data(i,:), 'LineWidth',2.0, 'Color', 'blue');
    plot(Time, dP_rot_wb_data(i,:), 'LineWidth',2.0, 'Color', 'cyan');
    plot(Time2, dP2_data(i,:), 'LineWidth',2.0, 'Color', 'magenta', 'LineStyle',':');
    ylabel('vel [$m/s$]', 'interpreter','latex', 'fontsize',15);
    axis tight;
    hold off;

    ax{3} = subplot(3,1,3);
    hold on;
    plot(Time, ddP_data(i,:), 'LineWidth',2.0, 'Color', 'blue');
    plot(Time, ddP_rot_wb_data(i,:), 'LineWidth',2.0, 'Color', 'cyan');
    plot(Time2, ddP2_data(i,:), 'LineWidth',2.0, 'Color', 'magenta', 'LineStyle',':');
    ylabel('accel [$m/s^2$]', 'interpreter','latex', 'fontsize',15);
    axis tight;
    hold off;

%     plotConstr(tau*x_pos_lim, pos_lb, pos_ub, tau*xeq_pos, pos_eq, ax{1}, i);
%     plotConstr(tau*x_vel_lim, vel_lb, vel_ub, tau*xeq_vel, vel_eq, ax{2}, i);
%     plotConstr(tau*x_accel_lim, accel_lb, accel_ub, tau*xeq_accel, accel_eq, ax{3}, i);
end

% ======================================================


function [Time, P_data, dP_data, ddP_data] = getGMPTrajectory(gmp, tau, y0, yg)

    Time = [];
    P_data = [];
    dP_data = [];
    ddP_data = [];

    p = y0;
    p_dot = zeros(size(p));
    p_ddot = zeros(size(p));
    
    gmp.setY0(y0);
    gmp.setGoal(yg);

    t = 0;
    dt = 0.002;

    while (true)

        x = t/tau;
        x_dot = 1/tau;
        
        if (x >= 1)
            x_dot = 0;
        end
        
        p_ref = gmp.getYd(x);
        p_ref_dot = gmp.getYdDot(x, x_dot);
        p_ref_ddot = gmp.getYdDDot(x, x_dot, 0);
        
        P = p_ref;
        p_dot = p_ref_dot;
        p_ddot = p_ref_ddot;

        Time = [Time t];
        P_data = [P_data p];
        dP_data = [dP_data p_dot];
        ddP_data = [ddP_data p_ddot];

        %p_ddot = p_ref_ddot + 30*(p_ref_dot - p_dot) + 100*(p_ref - p);

        t = t + dt;
        p = p + p_dot*dt;
        p_dot = p_dot + p_ddot*dt;

        if (x >= 1.0), break; end

    end

end



