clc;
% close all;
clear;

addpath('../../../../../matlab/lib/gmp_lib/');
import_gmp_lib();

addpath('../../../../../matlab/lib/io_lib/');
import_io_lib();

%% Load training data
fid = FileIO('data/constr_opt_pos_test_train_data.bin', FileIO.in);
Timed = fid.read('Timed');
Pd_data = fid.read('Pd_data');
dPd_data = fid.read('dPd_data');
ddPd_data = fid.read('ddPd_data');
fid.close();

Ts = Timed(2)-Timed(1);


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
% yg = ygd + [0.1; -0.1; 0.2]; view_ = [171.5301, -2.3630];
yg = ygd + [0.7; -0.7; 0.05];  view_ = [171.9421, -3.0690];


%% ======== Limits ==========
%            lower limit     upper limit
pos_lim = [[-1.2 -1.2 0.2]' [1.2 1.2 0.6]'];
vel_lim = [-0.3*ones(3,1) 0.3*ones(3,1)];  % lower and upper limit, same for all DoFs
accel_lim = [-0.4*ones(3,1) 0.4*ones(3,1)];

data = {};

% --------- Proportional scaling -----------
gmp.setScaleMethod(TrajScale_Prop(3));
[Time, P_data, dP_data, ddP_data] = getGMPTrajectory(gmp, tau, y0, yg);
data{length(data)+1} = ...
    struct('Time',Time, 'Pos',P_data, 'Vel',dP_data, 'Accel',ddP_data, 'linestyle',':', ...
        'color','blue', 'legend','prop', 'plot3D',true, 'plot2D',true);

% --------- Offline GMP-weights:VEL optimization -----------
[Time, P_data, dP_data, ddP_data] = offlineGMPweightsOpt(gmp, tau, y0, yg, pos_lim, vel_lim, accel_lim, false, true);
data{length(data)+1} = ...
    struct('Time',Time, 'Pos',P_data, 'Vel',dP_data, 'Accel',ddP_data, 'linestyle','-', ...
        'color',[0.85 0.33 0.1], 'legend','opt-w:vel', 'plot3D',true, 'plot2D',true);

% ---------- Offline GMP-trajectory optimization ------------
opt_pos = false;
opt_vel = true;
use_matlab_solver = 1;
[Time, P_data, dP_data, ddP_data] = offlineGMPtrajOpt(gmp, tau, y0, yg, pos_lim, vel_lim, accel_lim, opt_pos, opt_vel, use_matlab_solver);
data{length(data)+1} = ...
    struct('Time',Time, 'Pos',P_data, 'Vel',dP_data, 'Accel',ddP_data, 'linestyle',':', ...
        'color','green', 'legend','opt-y', 'plot3D',true, 'plot2D',true);
    
% ---------- Online GMP-trajectory optimization ------------
opt_pos = false;
opt_vel = true;
use_matlab_solver = 1;
[Time, P_data, dP_data, ddP_data] = onlineGMPtrajOpt(gmp, tau, y0, yg, pos_lim, vel_lim, accel_lim, opt_pos, opt_vel, use_matlab_solver);
data{length(data)+1} = ...
    struct('Time',Time, 'Pos',P_data, 'Vel',dP_data, 'Accel',ddP_data, 'linestyle',':', ...
    'color',[1, 0.41, 0.16], 'legend','onlineOpt-pos', 'plot3D',true, 'plot2D',true);


    
%% ======== Plot Results ==========

% plot 3D path
fig = figure;
fig.Position(3:4) = [815 716];
ax = axes();
hold on;
plot3(y0(1), y0(2), y0(3), 'LineWidth', 4, 'LineStyle','none', 'Color','green','Marker','o', 'MarkerSize',12);
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

title_ = {'x coordinate', 'y coordinate', 'z coordinate'};
label_font = 17;
ax_fontsize = 14;

for i=1:n_dof
    fig = figure;
    fig.Position(3:4) = [842 1110];

    ax = subplot(3,1,1);
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
    plot(0, y0(i), 'LineWidth', 4, 'LineStyle','none', 'Color','green','Marker','o', 'MarkerSize',10);
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
    hold on;
    for k=1:length(data)
        if (~data{k}.plot2D), continue; end
        plot(data{k}.Time, data{k}.Accel(i,:), 'LineWidth',2.5, 'LineStyle',data{k}.linestyle, 'Color',data{k}.color);
    end
    axis tight;
    % plot bounds
    plot(ax.XLim, [accel_lim(i,1) accel_lim(i,1)], 'LineWidth',2, 'LineStyle','--', 'Color',[1 0 1]);
    plot(ax.XLim, [accel_lim(i,2) accel_lim(i,2)], 'LineWidth',2, 'LineStyle','--', 'Color',[1 0 1]);
    ax.YLim = [ max(ax.YLim(1), 1.6*accel_lim(i,1)) min(ax.YLim(2), 1.6*accel_lim(i,2)) ];
    ylabel('accel [$m/s^2$]', 'interpreter','latex', 'fontsize',label_font);
    xlabel('time [$s$]', 'interpreter','latex', 'fontsize',label_font);
    ax.FontSize = ax_fontsize;
    hold off;

end



% ======================================================

function [Time, P_data, dP_data, ddP_data] = getOnlineOptGMPTrajectory2(gmp, tau, y0, yg, pos_lim, vel_lim, accel_lim, opt_pos, opt_vel)
    
    gmp2 = gmp.deepCopy();
    
    gmp2.setScaleMethod(TrajScale_Prop(3));
    gmp2.setY0(y0);
    gmp2.setGoal(yg);

    gmp_opt = GMP_Opt(gmp2);
    gmp_opt.setOptions(opt_pos, opt_vel, false, 0.1, 1, 0.1);
    gmp_opt.setMotionDuration(tau);

    % gmp_opt.optimize(100);
%     tic
%     gmp_opt.optimize2(0:0.01:1);
%     elaps_t = toc;
%     fprintf('===> Optimization finished! Elaps time: %f ms\n',elaps_t*1000);
%     fprintf([gmp_opt.getExitMsg() '\n']);
    
    Time = [];
    P_data = [];
    dP_data = [];
    ddP_data = [];

    p = y0;
    p_dot = zeros(size(p));
    p_ddot = zeros(size(p));
    
    gmp2.setY0(y0);
    gmp2.setGoal(yg);

    t = 0;
    dt = 0.002;
    x = t/tau;
    x_dot = 1/tau;
    x_ddot = 0;
    
    pos_lim = 10 * pos_lim;
    vel_lim = 10 * vel_lim;
    accel_lim = 10 * accel_lim;
    
    n_horiz = 4;
    x_horiz = 0:dt:(n_horiz-1)*dt;
    pos_lb = repmat( pos_lim(:,1), 1,n_horiz);
    pos_ub = repmat( pos_lim(:,2), 1,n_horiz);
    vel_lb = repmat( vel_lim(:,1)*ones(3,1), 1,n_horiz);
    vel_ub = repmat( vel_lim(:,2)*ones(3,1), 1,n_horiz);
    accel_lb = repmat( accel_lim(:,1)*ones(3,1), 1,n_horiz);
    accel_ub = repmat( accel_lim(:,2)*ones(3,1), 1,n_horiz);
    
    x_prev = 0;
    p_prev = gmp2.getYd(x);
    dp_prev = gmp2.getYdDot(x, x_dot);
    ddp_prev = gmp2.getYdDDot(x, x_dot, x_ddot);

    while (true)

        x = t/tau;
        x_dot = 1/tau;
        
        t/tau
        
        if (x >= 1)
            x_dot = 0;
        end
           
        % constraints
        gmp_opt.clearConstr();
        gmp_opt.setPosConstr(x_horiz, pos_lb, pos_ub, x_prev, p_prev);
%         gmp_opt.setVelConstr(x_horiz, vel_lb, vel_ub, x_prev, dp_prev);
%         gmp_opt.setAccelConstr(x_horiz, accel_lb, accel_ub, x_prev, ddp_prev);
%         gmp_opt.setPosConstr(x_horiz, pos_lb, pos_ub, [], []);
%         gmp_opt.setVelConstr(x_horiz, vel_lb, vel_ub, [], []);
%         gmp_opt.setAccelConstr(x_horiz, accel_lb, accel_ub, [], []);
        
        % optimize
        success = gmp_opt.optimize3(x_horiz, gmp);
        if (~success), warning([gmp_opt.getExitMsg() '\n']); end
        
        p_ref = gmp2.getYd(x);
        p_ref_dot = gmp2.getYdDot(x, x_dot);
        p_ref_ddot = gmp2.getYdDDot(x, x_dot, 0);
        
        p = p_ref;
        p_dot = p_ref_dot;
        p_ddot = p_ref_ddot;

        Time = [Time t];
        P_data = [P_data p];
        dP_data = [dP_data p_dot];
        ddP_data = [ddP_data p_ddot];

        %p_ddot = p_ref_ddot + 30*(p_ref_dot - p_dot) + 100*(p_ref - p);

        t = t + dt;
%         p = p + p_dot*dt;
%         p_dot = p_dot + p_ddot*dt;
        
        x_horiz = x_horiz + dt/tau;
        
        x_prev = x;
        p_prev = p;
        dp_prev = p_dot;
        ddp_prev = p_ddot;

        if (x >= 1.0), break; end

    end
    
end

function plot3Dbounds(ax, bounds)
    
    x1 = bounds(1,1);
    x2 = bounds(1,2);
    y1 = bounds(2,1);
    y2 = bounds(2,2);
    z1 = bounds(3,1);
    z2 = bounds(3,2);

%     patch( [x1 x1 x1 x1] , [y1 y1 y2 y2], [z1 z2 z2 z1], 'red', 'FaceAlpha',0.05, 'LineStyle','none', 'Parent',ax, 'HandleVisibility','off');
%     patch( [x2 x2 x2 x2] , [y1 y1 y2 y2], [z1 z2 z2 z1], 'red', 'FaceAlpha',0.05, 'LineStyle','none', 'Parent',ax, 'HandleVisibility','off');
    
    patch( [x1 x1 x1 x1] , [y1 y1 y2 y2], [z1 z2 z2 z1], 'red', 'FaceAlpha',0.05, 'LineStyle','none', 'Parent',ax, 'HandleVisibility','off');
    patch( [x2 x2 x2 x2] , [y1 y1 y2 y2], [z1 z2 z2 z1], 'red', 'FaceAlpha',0.05, 'LineStyle','none', 'Parent',ax, 'HandleVisibility','off');

    patch( [x1 x1 x2 x2] , [y1 y2 y2 y1], [z1 z1 z1 z1], 'red', 'FaceAlpha',0.05, 'LineStyle','none', 'Parent',ax, 'HandleVisibility','off');
    patch( [x1 x1 x2 x2] , [y1 y2 y2 y1], [z2 z2 z2 z2], 'red', 'FaceAlpha',0.05, 'LineStyle','none', 'Parent',ax, 'HandleVisibility','off');

%     patch( [x1 x1 x2 x2] , [y1 y1 y1 y1], [z1 z2 z2 z1], 'red', 'FaceAlpha',0.05, 'LineStyle','none', 'Parent',ax, 'HandleVisibility','off');
%     patch( [x1 x1 x2 x2] , [y2 y2 y2 y2], [z1 z2 z2 z1], 'red', 'FaceAlpha',0.05, 'LineStyle','none', 'Parent',ax, 'HandleVisibility','off');

end
