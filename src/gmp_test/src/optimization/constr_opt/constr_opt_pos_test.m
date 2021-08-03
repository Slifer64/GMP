clc;
close all;
clear;

addpath('../../../../matlab/lib/gmp_lib/');
import_gmp_lib();

addpath('../../../../matlab/lib/io_lib/');
import_io_lib();

%% Load training data
% load('data/constr_opt_pos_test_train_data.mat','Timed','Pd_data','dPd_data','ddPd_data');

% data = FileIO('train_data_pos_constr.bin', FileIO.in).readAll();
% Timed = data.Timed;
% Pd_data = data.Pd_data;

% fid = FileIO('data/constr_opt_pos_test_train_data.bin', bitor(FileIO.out, FileIO.trunc) );
% fid.write('Timed',Timed);
% fid.write('Pd_data',Pd_data);
% fid.write('dPd_data',dPd_data);
% fid.write('ddPd_data',ddPd_data);
% fid.close();

fid = FileIO('data/constr_opt_pos_test_train_data.bin', FileIO.in);
Timed = fid.read('Timed');
Pd_data = fid.read('Pd_data');
dPd_data = fid.read('dPd_data');
ddPd_data = fid.read('ddPd_data');
fid.close();

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


% ks = (yg - y0)/(ygd - yd0);
% ks

% gmp.setY0(y0);
% gmp.setGoal(yg);

% [Timed, Pd_data, dPd_data, ddPd_data] = getGMPTrajectory(gmp, tau, y0, yg);
% save('train_data_pos_constr.mat','Timed','Pd_data','dPd_data','ddPd_data');
% return;


%% ======== Limits ==========
%            lower limit     upper limit
pos_lim = [[-1.2 -1.2 0.2]' [1.2 1.2 0.6]'];
vel_lim = [-0.3 0.3];  % lower and upper limit, same for all DoFs
accel_lim = [-0.4 0.4];


%% ======== Generate trajectories ==========

% --------- Proportional scaling -----------
gmp.setScaleMethod(TrajScale_Prop(3));
[Time, P_data, dP_data, ddP_data] = getGMPTrajectory(gmp, tau, y0, yg);
data{1} = struct('Time',Time, 'Pos',P_data, 'Vel',dP_data, 'Accel',ddP_data, 'linestyle',':', ...
    'color','blue', 'legend','prop', 'plot3D',true, 'plot2D',true);

% --------- Rotational scaling -----------
traj_sc = TrajScale_Rot_wb();
traj_sc.setWorkBenchNormal([0; 0; 1]);
gmp.setScaleMethod(traj_sc);
[Time, P_data, dP_data, ddP_data] = getGMPTrajectory(gmp, tau, y0, yg);
data{2} = struct('Time',Time, 'Pos',P_data, 'Vel',dP_data, 'Accel',ddP_data, 'linestyle',':', ...
    'color','cyan', 'legend','rot-wb', 'plot3D',true, 'plot2D',true);

% --------- Demo -----------
% data{3} = struct('Time',Timed, 'Pos',Pd_data, 'Vel',dPd_data, 'Accel',ddPd_data, 'linestyle','--', 'color','magenta', 'legend','demo');
data{3} = struct('Time',Timed/kt, 'Pos',Pd_data, 'Vel',dPd_data*kt, 'Accel',ddPd_data*kt^2, 'linestyle','--', ...
    'color',[0.7 0 0], 'legend','demo', 'plot3D',true, 'plot2D',false);

% --------- Optimized DMP -> VEL -----------
[Time, P_data, dP_data, ddP_data] = getOptGMPTrajectory(gmp, tau, y0, yg, pos_lim, vel_lim, accel_lim, false, true);
data{4} = struct('Time',Time, 'Pos',P_data, 'Vel',dP_data, 'Accel',ddP_data, 'linestyle','-', ...
    'color','green', 'legend','opt-vel', 'plot3D',true, 'plot2D',true);

% --------- Optimized DMP -> POS -----------
[Time, P_data, dP_data, ddP_data] = getOptGMPTrajectory(gmp, tau, y0, yg, pos_lim, vel_lim, accel_lim, true, false);
data{5} = struct('Time',Time, 'Pos',P_data, 'Vel',dP_data, 'Accel',ddP_data, 'linestyle','-', ...
    'color',[0.85 0.33 0.1], 'legend','opt-pos', 'plot3D',true, 'plot2D',true);

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
    plot(ax.XLim, [vel_lim(1) vel_lim(1)], 'LineWidth',2, 'LineStyle','--', 'Color',[1 0 1]);
    plot(ax.XLim, [vel_lim(2) vel_lim(2)], 'LineWidth',2, 'LineStyle','--', 'Color',[1 0 1]);
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
    plot(ax.XLim, [accel_lim(1) accel_lim(1)], 'LineWidth',2, 'LineStyle','--', 'Color',[1 0 1]);
    plot(ax.XLim, [accel_lim(2) accel_lim(2)], 'LineWidth',2, 'LineStyle','--', 'Color',[1 0 1]);
    ylabel('accel [$m/s^2$]', 'interpreter','latex', 'fontsize',label_font);
    xlabel('time [$s$]', 'interpreter','latex', 'fontsize',label_font);
    ax.FontSize = ax_fontsize;
    hold off;

end

% ======================================================

function [Time, P_data, dP_data, ddP_data] = getOptGMPTrajectory(gmp, tau, y0, yg, pos_lim, vel_lim, accel_lim, opt_pos, opt_vel)
    
    gmp2 = gmp.deepCopy();
    
    gmp2.setScaleMethod(TrajScale_Prop(3));
    gmp2.setY0(y0);
    gmp2.setGoal(yg);

    gmp_opt = GMP_Opt(gmp2);
    gmp_opt.setOptions(opt_pos, opt_vel, false, 0.1, 1, 0.1);
    gmp_opt.setMotionDuration(tau);

    n_points = 100;
    % position constr
    gmp_opt.setPosBounds(pos_lim(:,1), pos_lim(:,2), n_points);
    gmp_opt.setPosConstr([],[],[], [0 1], [y0 yg]);
    % velocity constr
    gmp_opt.setVelBounds(vel_lim(1), vel_lim(2), n_points);
    gmp_opt.setVelConstr([], [], [], [0 1], zeros(3,2));
    % accel constr
    gmp_opt.setAccelBounds(accel_lim(1), accel_lim(2), n_points);
    gmp_opt.setAccelConstr([], [], [], [0 1], zeros(3,2));

    % gmp_opt.optimize(100);
    tic
    gmp_opt.optimize2(0:0.01:1);
    elaps_t = toc;
    fprintf('===> Optimization finished! Elaps time: %f ms\n',elaps_t*1000);
    fprintf([gmp_opt.getExitMsg() '\n']);
    
    [Time, P_data, dP_data, ddP_data] = getGMPTrajectory(gmp2, tau, y0, yg);
    
end

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


