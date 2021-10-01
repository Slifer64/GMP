clc;
% close all;
clear;

addpath('../../../../../matlab/lib/gmp_lib/');
import_gmp_lib();

addpath('../../../../../matlab/lib/io_lib/');
import_io_lib();

%% Load training data
fid = FileIO('data/train_data.bin', FileIO.in);
Timed = fid.read('Timed');
Pd_data = fid.read('Pd_data');
dPd_data = fid.read('dPd_data');
ddPd_data = fid.read('ddPd_data');
fid.close();

Ts = Timed(2)-Timed(1);

ind = [1 2 3];
Pd_data = Pd_data(ind,:);
dPd_data = dPd_data(ind,:);
ddPd_data = ddPd_data(ind,:);


%% initialize and train GMP
train_method = 'LS';
N_kernels = 40;
kernels_std_scaling = 1.5;
n_dof = size(Pd_data,1);
gmp = GMP(n_dof, N_kernels, kernels_std_scaling);
tic
offline_train_mse = gmp.train(train_method, Timed/Timed(end), Pd_data);
offline_train_mse
toc

taud = Timed(end);
yd0 = gmp.getYd(0); %Pd_data(1);
ygd = gmp.getYd(1); %Pd_data(end);

kt = 1.5; % temporal scaling
ks = diag([1 1 1]); % spatial scaling
tau = taud/kt;
y0 = yd0 + 0;
% yg = ks*(ygd - yd0) + y0;

%% case 1
% y_offset = [0.1; -0.1; 0.2];
% yg = ygd + y_offset(ind); view_ = [171.5301, -2.3630];

%% case 2
y_offset = [0.7; -0.7; 0.05];
yg = ygd + y_offset(ind);  view_ = [171.9421, -3.0690];


%% ======== Limits ==========
%            lower limit     upper limit
pos_lim = [[-1.2 -1.2 0.2]' [1.2 1.2 0.6]'];
pos_lim = pos_lim(ind,:);
vel_lim = 1*[-0.3*ones(n_dof,1) 0.3*ones(n_dof,1)];  % lower and upper limit, same for all DoFs
accel_lim = 3*[-0.4*ones(n_dof,1) 0.4*ones(n_dof,1)];

data = {};

% --------- Proportional scaling -----------
gmp.setScaleMethod(TrajScale_Prop(n_dof));
[Time, P_data, dP_data, ddP_data] = getGMPTrajectory(gmp, tau, y0, yg);
data{length(data)+1} = ...
    struct('Time',Time, 'Pos',P_data, 'Vel',dP_data, 'Accel',ddP_data, 'linestyle',':', ...
        'color','blue', 'legend','prop', 'plot3D',true, 'plot2D',true);

opt_pos = 1;
opt_vel = 0;
qp_solver_type = 1; % matlab-quadprog:0 , osqp:1, Goldfarb-Idnani: 2
opt_type = 'pos';
if (opt_vel), opt_type = 'vel'; end
% 
% % --------- Offline GMP-weights optimization -----------
% [Time, P_data, dP_data, ddP_data] = offlineGMPweightsOpt(gmp, tau, y0, yg, pos_lim, vel_lim, accel_lim, opt_pos, opt_vel);
% data{length(data)+1} = ...
%     struct('Time',Time, 'Pos',P_data, 'Vel',dP_data, 'Accel',ddP_data, 'linestyle',':', ...
%         'color',[0.64,0.08,0.18], 'legend',['opt-w:' opt_type], 'plot3D',true, 'plot2D',true);
% 
% % ---------- Online GMP-weights optimization ------------
% [Time, P_data, dP_data, ddP_data] = onlineGMPweightsOpt(gmp, tau, y0, yg, pos_lim, 1*vel_lim, accel_lim, opt_pos, opt_vel, qp_solver_type);
% data{length(data)+1} = ...
%     struct('Time',Time, 'Pos',P_data, 'Vel',dP_data, 'Accel',ddP_data, 'linestyle','-', ...
%     'color',[1, 0.41, 0.16], 'legend',['opt-w:' opt_type '(online)'], 'plot3D',true, 'plot2D',true);

% % ---------- GMP-MPC optimization ------------
% [Time, P_data, dP_data, ddP_data] = gmpMpcOpt(gmp, tau, y0, yg, pos_lim, 1*vel_lim, accel_lim, opt_pos, opt_vel, qp_solver_type);
% data{length(data)+1} = ...
%     struct('Time',Time, 'Pos',P_data, 'Vel',dP_data, 'Accel',ddP_data, 'linestyle','-', ...
%     'color',[0.72 0.27 1], 'legend',['gmp-mpc:' opt_type], 'plot3D',true, 'plot2D',true);

[Time, P_data, dP_data, ddP_data] = gmpMpcOptCpp(gmp, tau, y0, yg, pos_lim, 1*vel_lim, accel_lim, 'out');
data{length(data)+1} = ...
    struct('Time',Time, 'Pos',P_data, 'Vel',dP_data, 'Accel',ddP_data, 'linestyle','-', ...
    'color',[1 0.07 0.65], 'legend',['gmp-mpc:' opt_type], 'plot3D',true, 'plot2D',true);


% % ---------- Offline GMP-trajectory optimization ------------
% [Time, P_data, dP_data, ddP_data] = offlineGMPtrajOpt(gmp, tau, y0, yg, pos_lim, vel_lim, accel_lim, opt_pos, opt_vel, qp_solver_type);
% data{length(data)+1} = ...
%     struct('Time',Time, 'Pos',P_data, 'Vel',dP_data, 'Accel',ddP_data, 'linestyle',':', ...
%         'color','green', 'legend',['opt-traj:' opt_type], 'plot3D',true, 'plot2D',true);
%      
% % ---------- Online GMP-trajectory optimization ------------
% [Time, P_data, dP_data, ddP_data] = onlineGMPtrajOpt(gmp, tau, y0, yg, pos_lim + 0*[-1 1], 1*vel_lim, 1*accel_lim, opt_pos, opt_vel, qp_solver_type);
% data{length(data)+1} = ...
%     struct('Time',Time, 'Pos',P_data, 'Vel',dP_data, 'Accel',ddP_data, 'linestyle','-', ...
%     'color',[0 0.7 0], 'legend',['opt-traj:' opt_type '(online)'], 'plot3D',true, 'plot2D',true);

% % ---------- GMP with repulsive forces ------------
% [Time, P_data, dP_data, ddP_data] = gmpWithRepulsiveForces(gmp, tau, y0, yg, pos_lim, vel_lim, accel_lim);
% data{length(data)+1} = ...
%     struct('Time',Time, 'Pos',P_data, 'Vel',dP_data, 'Accel',ddP_data, 'linestyle','-', ...
%     'color',[0.93 0.69 0.13], 'legend',['gmp-repforce'], 'plot3D',true, 'plot2D',true);

% dt = Time(2)-Time(1);
% dP1_data = [diff(P_data) 0] / dt;
% ddP1_data = [diff(dP1_data) 0] / dt;
% ddP2_data = [diff(dP_data) 0] / dt;
% 
% figure;
% ax1 = subplot(3,1,1); hold on;
% plot(Time, P_data, 'LineWidth',2, 'LineStyle','-');
% ax2 = subplot(3,1,2); hold on;
% plot(Time, dP_data, 'LineWidth',2, 'LineStyle','-', 'Color','blue');
% plot(Time, dP1_data, 'LineWidth',2, 'LineStyle','--', 'Color','magenta');
% ax3 = subplot(3,1,3); hold on;
% plot(Time, ddP_data, 'LineWidth',2, 'LineStyle','-', 'Color','blue');
% plot(Time, ddP1_data, 'LineWidth',2, 'LineStyle','--', 'Color','magenta');
% plot(Time, ddP2_data, 'LineWidth',2, 'LineStyle',':', 'Color','green');
% linkaxes([ax1 ax2 ax3],'x');    

%% ======== Plot Results ==========

% plot 3D path
if (n_dof == 3)
    fig = figure;
    fig.Position(3:4) = [815 716];
    ax = axes();
    hold on;
    plot3(y0(1), y0(2), y0(3), 'LineWidth', 4, 'LineStyle','none', 'Color',[0.47,0.67,0.19],'Marker','o', 'MarkerSize',12);
    plot3(yg(1), yg(2), yg(3), 'LineWidth', 4, 'LineStyle','none', 'Color','red','Marker','x', 'MarkerSize',12);
    plot3(ygd(1), ygd(2), ygd(3), 'LineWidth', 4, 'LineStyle','none', 'Color','magenta','Marker','x', 'MarkerSize',12);
    legend_ = {};
    for k=1:length(data)
        if (~data{k}.plot3D), continue; end
        plot3(data{k}.Pos(1,:), data{k}.Pos(2,:), data{k}.Pos(3,:), 'LineWidth', 3, 'LineStyle',data{k}.linestyle, 'Color',data{k}.color);
        legend_ = [legend_ data{k}.legend];
    end
    plot3([ygd(1) yg(1)], [ygd(2) yg(2)], [ygd(3) yg(3)], 'LineWidth', 1, 'LineStyle','--', 'Color',[1 0 1 0.5]);
    legend(['$p_0$','$g$','$g_d$' legend_], 'interpreter','latex', 'fontsize',17, 'Position',[0.8174 0.6641 0.1531 0.3144]);
    axis tight;
    x_lim = ax.XLim + 0.05*[-1 1];
    y_lim = ax.YLim + 0.05*[-1 1];
    z_lim = ax.ZLim + 0.05*[-1 1];
    plot3Dbounds(ax, pos_lim)
    view(view_);
    grid on;
    ax.XLim=x_lim; ax.YLim=y_lim; ax.ZLim=z_lim;
    ax.FontSize = 14;
    hold off;
end

title_ = {'x coordinate', 'y coordinate', 'z coordinate'};
label_font = 17;
ax_fontsize = 14;

for i=1:n_dof
    fig = figure;
    fig.Position(3:4) = [842 1110];

    ax = subplot(3,1,1);
    ax_vec = [ax];
    hold on;
    % plot position trajectory
    legend_ = {};
    for k=1:length(data)
        if (~data{k}.plot2D), continue; end
        plot(data{k}.Time, data{k}.Pos(i,:), 'LineWidth',2.5, 'LineStyle',data{k}.linestyle, 'Color',data{k}.color);
        legend_ = [legend_ data{k}.legend];
    end
    axis tight;
    % plot start and final positions
    plot(0, y0(i), 'LineWidth', 4, 'LineStyle','none', 'Color',[0.47,0.67,0.19],'Marker','o', 'MarkerSize',10);
    plot(tau, yg(i), 'LineWidth', 4, 'LineStyle','none', 'Color','red','Marker','x', 'MarkerSize',10);
    plot(tau, ygd(i), 'LineWidth', 4, 'LineStyle','none', 'Color','magenta','Marker','x', 'MarkerSize',10); 
    % plot bounds
    plot(ax.XLim, [pos_lim(i,1) pos_lim(i,1)], 'LineWidth',2, 'LineStyle','--', 'Color',[1 0 1]);
    plot(ax.XLim, [pos_lim(i,2) pos_lim(i,2)], 'LineWidth',2, 'LineStyle','--', 'Color',[1 0 1]);
    % labels, title ...
    ylabel('pos [$m$]', 'interpreter','latex', 'fontsize',label_font);
%     title(title_{i}, 'interpreter','latex', 'fontsize',18);
    legend(legend_, 'interpreter','latex', 'fontsize',17, 'Position',[0.2330 0.9345 0.5520 0.0294], 'Orientation', 'horizontal');
    ax.FontSize = ax_fontsize;
    hold off;

    ax = subplot(3,1,2);
    ax_vec = [ax_vec ax];
    hold on;
    for k=1:length(data) 
        if (~data{k}.plot2D), continue; end
        plot(data{k}.Time, data{k}.Vel(i,:), 'LineWidth',2.5, 'LineStyle',data{k}.linestyle, 'Color',data{k}.color);
    end
    axis tight;
    % plot bounds
    plot(ax.XLim, [vel_lim(i,1) vel_lim(i,1)], 'LineWidth',2, 'LineStyle','--', 'Color',[1 0 1]);
    plot(ax.XLim, [vel_lim(i,2) vel_lim(i,2)], 'LineWidth',2, 'LineStyle','--', 'Color',[1 0 1]);
    ylabel('vel [$m/s$]', 'interpreter','latex', 'fontsize',label_font);
    ax.FontSize = ax_fontsize;
    hold off;

    ax = subplot(3,1,3);
    ax_vec = [ax_vec ax];
    hold on;
    for k=1:length(data)
        if (~data{k}.plot2D), continue; end
        plot(data{k}.Time, data{k}.Accel(i,:), 'LineWidth',2.5, 'LineStyle',data{k}.linestyle, 'Color',data{k}.color);
    end
    axis tight;
    % plot bounds
    plot(ax.XLim, [accel_lim(i,1) accel_lim(i,1)], 'LineWidth',2, 'LineStyle','--', 'Color',[1 0 1]);
    plot(ax.XLim, [accel_lim(i,2) accel_lim(i,2)], 'LineWidth',2, 'LineStyle','--', 'Color',[1 0 1]);
    %ax.YLim = [ max(ax.YLim(1), 8*accel_lim(i,1)) min(ax.YLim(2), 8*accel_lim(i,2)) ];
    ylabel('accel [$m/s^2$]', 'interpreter','latex', 'fontsize',label_font);
    xlabel('time [$s$]', 'interpreter','latex', 'fontsize',label_font);
    ax.FontSize = ax_fontsize;
    hold off;
    
    linkaxes(ax_vec,'x');

end

% ======================================================

function plot3Dbounds(ax, bounds)
    
    x1 = bounds(1,1);
    x2 = bounds(1,2);
    y1 = bounds(2,1);
    y2 = bounds(2,2);
    z1 = bounds(3,1);
    z2 = bounds(3,2);

%     patch( [x1 x1 x1 x1] , [y1 y1 y2 y2], [z1 z2 z2 z1], 'red', 'FaceAlpha',0.05, 'LineStyle','none', 'Parent',ax, 'HandleVisibility','off');
%     patch( [x2 x2 x2 x2] , [y1 y1 y2 y2], [z1 z2 z2 z1], 'red', 'FaceAlpha',0.05, 'LineStyle','none', 'Parent',ax, 'HandleVisibility','off');
    
    X = [x1 x1 x1 x1; x2 x2 x2 x2; x1 x1 x2 x2; x1 x1 x2 x2]';
    Y = [y1 y1 y2 y2; y1 y1 y2 y2; y1 y2 y2 y1; y1 y2 y2 y1]';
    Z = [z1 z2 z2 z1; z1 z2 z2 z1; z1 z1 z1 z1; z2 z2 z2 z2]';
    
    patch( X , Y, Z, 'red', 'FaceAlpha',0.05, 'LineStyle','none', 'Parent',ax, 'HandleVisibility','off');

%     patch( [x1 x1 x2 x2] , [y1 y1 y1 y1], [z1 z2 z2 z1], 'red', 'FaceAlpha',0.05, 'LineStyle','none', 'Parent',ax, 'HandleVisibility','off');
%     patch( [x1 x1 x2 x2] , [y2 y2 y2 y2], [z1 z2 z2 z1], 'red', 'FaceAlpha',0.05, 'LineStyle','none', 'Parent',ax, 'HandleVisibility','off');

end
