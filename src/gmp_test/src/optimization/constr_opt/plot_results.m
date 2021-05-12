clc;
close all;
clear;

addpath('../../../../matlab/lib/io_lib/');
import_io_lib();

%% Load results
fid = FileIO('results.bin', FileIO.in);
Time = fid.read('Time');
P_data = fid.read('P_data');
dP_data = fid.read('dP_data');
ddP_data = fid.read('ddP_data');
% -----------------------------
Time2 = fid.read('Time2');
P_data2 = fid.read('P_data2');
dP_data2 = fid.read('dP_data2');
ddP_data2 = fid.read('ddP_data2');
% -----------------------------
fid.close();

%% Load constraints
fid = FileIO('constraints.bin', FileIO.in);
t_pos_lim = fid.read('t_pos_lim');
pos_lb = fid.read('pos_lb');
pos_ub = fid.read('pos_ub');
teq_pos = fid.read('teq_pos');
pos_eq = fid.read('pos_eq');
t_vel_lim = fid.read('t_vel_lim');
vel_lb = fid.read('vel_lb');
vel_ub = fid.read('vel_ub');
teq_vel = fid.read('teq_vel');
vel_eq = fid.read('vel_eq');
t_accel_lim = fid.read('t_accel_lim');
accel_lb = fid.read('accel_lb');
accel_ub = fid.read('accel_ub');
teq_accel = fid.read('teq_accel');
accel_eq = fid.read('accel_eq');
fid.close();


%% plot results

% plot 3D path
figure;
hold on;
plot3(P_data(1,:), P_data(2,:), P_data(3,:), 'LineWidth', 2, 'LineStyle','-', 'Color','blue');
plot3(P_data2(1,:), P_data2(2,:), P_data2(3,:), 'LineWidth', 2, 'LineStyle',':', 'Color',[0.85 0.33 0.1]);
legend({'$gmp$', '$k_s * demo$'}, 'interpreter','latex', 'fontsize',17);
hold off;

% --------------------------

n_dof = size(P_data,1);

for i=1:n_dof
    figure;
    ax = cell(3,1);

    ax{1} = subplot(3,1,1);
    hold on;
    plot(Time, P_data(i,:), 'LineWidth',2.0, 'Color', 'blue');
    plot(Time2, P_data2(i,:), 'LineWidth',2.0, 'Color', 'magenta', 'LineStyle',':');
    ylabel('pos [$m$]', 'interpreter','latex', 'fontsize',15);
    %title(['temporal scale: $' num2str(kt) '$     ,     spatial scale: $' num2str(det(ks)) '$'], 'interpreter','latex', 'fontsize',18);
    scatter(nan, nan, 'MarkerEdgeColor',[0 0.7 0], 'Marker','o', 'MarkerEdgeAlpha',0.6, 'LineWidth',4, 'SizeData', 100); % dummy plot for legend
    scatter(nan, nan, 'MarkerEdgeColor',[1 0 0], 'Marker','^', 'MarkerEdgeAlpha',0.3, 'LineWidth',2, 'SizeData', 100); % dummy plot for legend
    scatter(nan, nan, 'MarkerEdgeColor',[1 0 0], 'Marker','v', 'MarkerEdgeAlpha',0.3, 'LineWidth',2, 'SizeData', 100); % dummy plot for legend
    % legend({'GMP','ref','eq-constr','low-bound','upper-bound'}, 'interpreter','latex', 'fontsize',15);
    axis tight;
    hold off;

    ax{2} = subplot(3,1,2);
    hold on;
    plot(Time, dP_data(i,:), 'LineWidth',2.0, 'Color', 'blue');
    plot(Time2, dP_data2(i,:), 'LineWidth',2.0, 'Color', 'magenta', 'LineStyle',':');
    ylabel('vel [$m/s$]', 'interpreter','latex', 'fontsize',15);
    axis tight;
    hold off;

    ax{3} = subplot(3,1,3);
    hold on;
    plot(Time, ddP_data(i,:), 'LineWidth',2.0, 'Color', 'blue');
    plot(Time2, ddP_data2(i,:), 'LineWidth',2.0, 'Color', 'magenta', 'LineStyle',':');
    ylabel('accel [$m/s^2$]', 'interpreter','latex', 'fontsize',15);
    axis tight;
    hold off;

    plotConstr(t_pos_lim, pos_lb, pos_ub, teq_pos, pos_eq, ax{1}, i);
    plotConstr(t_vel_lim, vel_lb, vel_ub, teq_vel, vel_eq, ax{2}, i);
    plotConstr(t_accel_lim, accel_lb, accel_ub, teq_accel, accel_eq, ax{3}, i);
end

%% =======================================================
%% =======================================================

function plotConstr(t, lb, ub, teq, feq, ax, i_dof)

    hold(ax, 'on');
    
    if (~isempty(teq))
        eq_sc = scatter(teq, feq(i_dof,:), 'MarkerEdgeColor',[0 0.7 0], 'Marker','*', ...
                'MarkerEdgeAlpha',0.6, 'LineWidth',2, 'SizeData', 100, 'Parent',ax);
    end
    
    if (~isempty(t))
        less_sc = scatter(t, lb(i_dof,:), 'MarkerEdgeColor',[1 0 0], 'Marker','.', ...
                'MarkerEdgeAlpha',0.3, 'LineWidth',1, 'SizeData', 100, 'Parent',ax);
        %less_sc = plot(t, lb(i_dof,:), 'Color',[1 0 0 0.7], 'LineStyle','--', 'LineWidth',2, 'Parent',ax);

        gr_sc = scatter(t, ub(i_dof,:), 'MarkerEdgeColor',[1 0 0], 'Marker','.', ...
                'MarkerEdgeAlpha',0.3, 'LineWidth',1, 'SizeData', 100, 'Parent',ax);
%         gr_sc = plot(t, ub(i_dof,:), 'Color',[0 0.8 0 0.7], 'LineStyle','--', 'LineWidth',2, 'Parent',ax);
    end

end
