clc;
close all;
clear;

% pos_lim = [[-1.2 -1.2 0.1]' [1.2 1.2 0.6]'];
% x1 = pos_lim(1,1);
% x2 = pos_lim(1,2);
% y1 = pos_lim(2,1);
% y2 = pos_lim(2,2);
% z1 = pos_lim(3,1);
% z2 = pos_lim(3,2);
% % 
% % X = [x1 x1 x1 x1 x2 x2 x2 x2];
% % Y = [y1 y2 y2 y1]
% % Z = [z1 z1 z2 z2]
% 
% figure;
% hold on;
% patch( [x1 x1 x1 x1] , [y1 y1 y2 y2], [z1 z2 z2 z1], 'red', 'FaceAlpha',0.05, 'LineStyle','none');
% patch( [x2 x2 x2 x2] , [y1 y1 y2 y2], [z1 z2 z2 z1], 'red', 'FaceAlpha',0.05, 'LineStyle','none');
% 
% patch( [x1 x1 x2 x2] , [y1 y2 y2 y1], [z1 z1 z1 z1], 'red', 'FaceAlpha',0.05, 'LineStyle','none');
% patch( [x1 x1 x2 x2] , [y1 y2 y2 y1], [z2 z2 z2 z2], 'red', 'FaceAlpha',0.05, 'LineStyle','none');
%
% patch( [x1 x1 x2 x2] , [y1 y1 y1 y1], [z1 z2 z2 z1], 'red', 'FaceAlpha',0.05, 'LineStyle','none');
% patch( [x1 x1 x2 x2] , [y2 y2 y2 y2], [z1 z2 z2 z1], 'red', 'FaceAlpha',0.05, 'LineStyle','none');
% 
% view(9.6, 19.2);
% return


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
pos_lim = [[-1.2 -1.2 0.2]' [1.2 1.2 0.6]'];
vel_lim = [-0.3 0.3];
accel_lim = [-0.4 0.4];


%% ======== Generate trajectories ==========

% --------- Proportional scaling -----------
gmp.setScaleMethod(TrajScale_Prop(3));
[Time, P_data, dP_data, ddP_data] = getGMPTrajectory(gmp, tau, y0, yg);
data{1} = struct('Time',Time, 'Pos',P_data, 'Vel',dP_data, 'Accel',ddP_data, 'linestyle','-', 'color','blue', 'legend','prop');

% --------- Rotational scaling -----------
gmp.setScaleMethod(TrajScale_Rot_wb());
[Time, P_data, dP_data, ddP_data] = getGMPTrajectory(gmp, tau, y0, yg);
data{2} = struct('Time',Time, 'Pos',P_data, 'Vel',dP_data, 'Accel',ddP_data, 'linestyle','-', 'color','cyan', 'legend','rot-wb');

% --------- Demo -----------
data{3} = struct('Time',Timed/kt, 'Pos',Pd_data, 'Vel',dPd_data*kt, 'Accel',ddPd_data*kt^2, 'linestyle','--', 'color','magenta', 'legend','demo');

% --------- Optimized DMP -> VEL -----------
[Time, P_data, dP_data, ddP_data] = getOptGMPTrajectory(gmp, tau, y0, yg, pos_lim, vel_lim, accel_lim, false, true);
data{4} = struct('Time',Time, 'Pos',P_data, 'Vel',dP_data, 'Accel',ddP_data, 'linestyle','-', 'color','green', 'legend','opt-vel');

[Time, P_data, dP_data, ddP_data] = getOptGMPTrajectory(gmp, tau, y0, yg, pos_lim, vel_lim, accel_lim, true, false);
data{5} = struct('Time',Time, 'Pos',P_data, 'Vel',dP_data, 'Accel',ddP_data, 'linestyle','-', 'color',[0 0.5 0], 'legend','opt-pos');


%% ======== Plot Results ==========

% plot 3D path
fig = figure;
fig.Position(3:4) = [815 716];
ax = axes();
hold on;
plot3(y0(1), y0(2), y0(3), 'LineWidth', 4, 'LineStyle','none', 'Color','green','Marker','o', 'MarkerSize',8);
plot3(yg(1), yg(2), yg(3), 'LineWidth', 4, 'LineStyle','none', 'Color','red','Marker','x', 'MarkerSize',8);
plot3(ygd(1), ygd(2), ygd(3), 'LineWidth', 4, 'LineStyle','none', 'Color','magenta','Marker','x', 'MarkerSize',8);
legend_ = {};
for k=1:length(data)
    dat = data{k};
    plot3(dat.Pos(1,:), dat.Pos(2,:), dat.Pos(3,:), 'LineWidth', 2, 'LineStyle',dat.linestyle, 'Color',dat.color);
    legend_ = [legend_ dat.legend];
end
plot3([ygd(1) yg(1)], [ygd(2) yg(2)], [ygd(3) yg(3)], 'LineWidth', 1, 'LineStyle','--', 'Color',[1 0 1 0.5]);
legend(['$p_0$','$g$','$g_d$' legend_], 'interpreter','latex', 'fontsize',17, 'Position',[0.8014 0.6278 0.1531 0.3144]);
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




return

for i=1:n_dof
    figure;
    ax = cell(3,1);

    ax{1} = subplot(3,1,1);
    hold on;
    for k=1:length(data) 
        plot(data{k}.Time, data{k}.Pos(i,:), 'LineWidth',2.0, 'LineStyle',data{k}.linestyle, 'Color',data{k}.color);
    end
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
    for k=1:length(data) 
        plot(data{k}.Time, data{k}.Vel(i,:), 'LineWidth',2.0, 'LineStyle',data{k}.linestyle, 'Color',data{k}.color);
    end
    ylabel('vel [$m/s$]', 'interpreter','latex', 'fontsize',15);
    axis tight;
    hold off;

    ax{3} = subplot(3,1,3);
    hold on;
    for k=1:length(data) 
        plot(data{k}.Time, data{k}.Accel(i,:), 'LineWidth',2.0, 'LineStyle',data{k}.linestyle, 'Color',data{k}.color);
    end
    ylabel('accel [$m/s^2$]', 'interpreter','latex', 'fontsize',15);
    axis tight;
    hold off;

%     plotConstr(tau*x_pos_lim, pos_lb, pos_ub, tau*xeq_pos, pos_eq, ax{1}, i);
%     plotConstr(tau*x_vel_lim, vel_lb, vel_ub, tau*xeq_vel, vel_eq, ax{2}, i);
%     plotConstr(tau*x_accel_lim, accel_lb, accel_ub, tau*xeq_accel, accel_eq, ax{3}, i);
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
    toc
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

    patch( [x1 x1 x2 x2] , [y1 y2 y2 y1], [z1 z1 z1 z1], 'red', 'FaceAlpha',0.05, 'LineStyle','none', 'Parent',ax, 'HandleVisibility','off');
    patch( [x1 x1 x2 x2] , [y1 y2 y2 y1], [z2 z2 z2 z2], 'red', 'FaceAlpha',0.05, 'LineStyle','none', 'Parent',ax, 'HandleVisibility','off');

%     patch( [x1 x1 x2 x2] , [y1 y1 y1 y1], [z1 z2 z2 z1], 'red', 'FaceAlpha',0.05, 'LineStyle','none', 'Parent',ax, 'HandleVisibility','off');
%     patch( [x1 x1 x2 x2] , [y2 y2 y2 y2], [z1 z2 z2 z1], 'red', 'FaceAlpha',0.05, 'LineStyle','none', 'Parent',ax, 'HandleVisibility','off');

end


