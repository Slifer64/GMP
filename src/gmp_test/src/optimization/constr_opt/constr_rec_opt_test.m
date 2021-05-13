clc;
close all;
clear;

addpath('../../../../matlab/lib/gmp_lib/');
import_gmp_lib();

addpath('../../../../matlab/lib/io_lib/');
import_io_lib();

%% Load training data
fid = FileIO('pos_data.bin', FileIO.in);
Timed = fid.read('Timed');
Pd_data = fid.read('Pd_data');
dPd_data = fid.read('dPd_data');
ddPd_data = fid.read('ddPd_data');
fid.close();

Ts = Timed(2)-Timed(1);

%% initialize and train GMP
train_method = 'LS';
N_kernels = 30;
kernels_std_scaling = 1.5;
n_dof = size(Pd_data,1);
gmp = GMP_nDoF(n_dof, N_kernels, kernels_std_scaling);
tic
offline_train_mse = gmp.train(train_method, Timed/Timed(end), Pd_data);
offline_train_mse
toc

% gmp.setScaleMethod(TrajScale.ROT_MIN_SCALE);

taud = Timed(end);
yd0 = gmp.getYd(0); %Pd_data(1);
ygd = gmp.getYd(1); %Pd_data(end);

kt = 3; % temporal scaling
ks = diag([1.3 1.4 1.5]); % spatial scaling
tau = taud/kt;
y0 = yd0 + 0.15;
yg = ks*(ygd - yd0) + y0;
gmp.setY0(y0);
gmp.setGoal(yg);


%% calculate scaled demo trajectory
T_sc = gmp.getScaling();
Time2 = Timed / kt;
P_data2 = T_sc*(Pd_data - Pd_data(:,1)) + y0;
dP_data2 = kt*T_sc*dPd_data;
ddP_data2 = (kt^2)*T_sc*ddPd_data;


%% set up constraints
n_data = length(Time2);
P_max = max([P_data2(:,1) P_data2(:,end)]')' + 0.1;
P_min = min([P_data2(:,1) P_data2(:,end)]')' - 0.1;
ind = 1:round(n_data/100):n_data;

x_lim = Time2(ind)/Time2(end);
Pos_up_lim = P_data2(:,ind) + 0.3;
Pos_low_lim = P_data2(:,ind) - 0.3;
for i=1:3
    Pos_up_lim(i, Pos_up_lim(i,:)>P_max(i)) = P_max(i);
    Pos_low_lim(i, Pos_low_lim(i,:)<P_min(i)) = P_min(i);
    
    ind2 = Pos_low_lim(i,:)> Pos_up_lim(i,:)-0.3;
    Pos_low_lim(i, ind2) = Pos_up_lim(i,ind2) - 0.3;
end

x_pos_lim = x_lim;
pos_lb = Pos_low_lim;
pos_ub = Pos_up_lim;

xeq_pos = [0 1]; % [];
pos_eq = [y0 yg];%[]; %

x_vel_lim = 0:0.01:1;
n_vel_constr = length(x_vel_lim);
vel_lb = -1*ones(n_dof,n_vel_constr);
vel_ub = 1*ones(n_dof,n_vel_constr);

xeq_vel = [0 1];%[]; %
vel_eq = zeros(3,2);%[]; %

x_accel_lim = 0:0.01:1;
n_accel_constr = length(x_accel_lim);
accel_lb = -2*ones(n_dof,n_accel_constr);
accel_ub = 2*ones(n_dof,n_accel_constr);

xeq_accel = [0 1];%[]; % 
accel_eq = zeros(3,2);%[]; % 

% %% Write constraints to file
% if (false) % set to 'true' to write constraints to file
%     fid = FileIO ('constraints.bin', bitor(FileIO.out, FileIO.trunc) );
%     % -----------------------------
%     fid.write('t_pos_lim', tau*x_pos_lim);
%     fid.write('pos_lb', pos_lb);
%     fid.write('pos_ub', pos_ub);
%     fid.write('teq_pos', tau*xeq_pos);
%     fid.write('pos_eq', pos_eq);
%     fid.write('t_vel_lim', tau*x_vel_lim);
%     fid.write('vel_lb', vel_lb);
%     fid.write('vel_ub', vel_ub);
%     fid.write('teq_vel', tau*xeq_vel);
%     fid.write('vel_eq', vel_eq);
%     fid.write('t_accel_lim', tau*x_accel_lim);
%     fid.write('accel_lb', accel_lb);
%     fid.write('accel_ub', accel_ub);
%     fid.write('teq_accel', tau*xeq_accel);
%     fid.write('accel_eq', accel_eq);
%     % -----------------------------
%     fid.close();
% end

%% Solve constr-optimization problem on-line
gmp_opt = GMP_nDoF_Opt(gmp);
gmp_opt.setOptions(true, true, false, 0.1, 1, 0.1);

Time = [];
P_data = [];
dP_data = [];
ddP_data = [];

Elaps_Time = [];

p = y0;
p_dot = 0;
p_ddot = 0;

t = 0;
dt = 0.01;

dx = 0.01;
x_win = 16 / sqrt(2*gmp.h(1)); % 4*sigma
x_hor = 0:dx:x_win;
n_hor = length(x_hor);

pos_lb_i = zeros(n_dof, n_hor);
pos_ub_i = zeros(n_dof, n_hor);

vel_lb_i = -1 * ones(n_dof, n_hor);
vel_ub_i = 1 * ones(n_dof, n_hor);
accel_lb_i = -2 * ones(n_dof, n_hor);
accel_ub_i = 2 * ones(n_dof, n_hor);

x = 0;
x_dot = 1 / tau;
x_ddot = 0;

W0 = gmp.W;

while (true)
    
    for k=1:n_dof
        pos_lb_i(k,:) = interp1(x_pos_lim, pos_lb(k,:), x_hor);
        pos_ub_i(k,:) = interp1(x_pos_lim, pos_ub(k,:), x_hor);
    end
    
    tic;
    
%     % set equality constraint to current trajectory point
%     xeq_pos_i = x_hor(1);
%     pos_eq_i = gmp.getYd(x);
%     xeq_vel_i = x_hor(1);
%     vel_eq_i = gmp.getYdDot(x_hor(1), x_dot);
%     xeq_accel_i = x_hor(1);
%     accel_eq_i = gmp.getYdDDot(x, x_dot, x_ddot);
    
    gmp.W = W0;
    gmp_opt.setMotionDuration(tau);
    gmp_opt.setPosConstr(x_hor, pos_lb_i, pos_ub_i, [], []);
    gmp_opt.setVelConstr(x_hor, vel_lb_i, vel_ub_i, [], []);
    gmp_opt.setAccelConstr(x_hor, accel_lb_i, accel_ub_i, [], []);
    success = gmp_opt.optimize2(x_hor);

    if (~success), warning('Opt failed...'); end
    
    x = t/tau;
    p = gmp.getYd(x);
    p_dot = gmp.getYdDot(x, x_dot);
    p_ddot = gmp.getYdDDot(x, x_dot, 0);
    
    Elaps_Time = [Elaps_Time toc];
    
    [t Elaps_Time(end)]

    Time = [Time t];
    P_data = [P_data p];
    dP_data = [dP_data p_dot];
    ddP_data = [ddP_data p_ddot];

    t = t + dt;
    
    x_hor = x_hor + x_dot*dt;

    if (t >= 2), break; end

end


%% plot results

% plot 3D path
figure;
hold on;
plot3(P_data(1,:), P_data(2,:), P_data(3,:), 'LineWidth', 2, 'LineStyle','-', 'Color','blue');
plot3(P_data2(1,:), P_data2(2,:), P_data2(3,:), 'LineWidth', 2, 'LineStyle',':', 'Color',[0.85 0.33 0.1]);
legend({'$gmp$', '$k_s * demo$'}, 'interpreter','latex', 'fontsize',17);
hold off;

% --------------------------

for i=1:n_dof
    figure;
    ax = cell(3,1);

    ax{1} = subplot(3,1,1);
    hold on;
    plot(Time, P_data(i,:), 'LineWidth',2.0, 'Color', 'blue');
    plot(Time2, P_data2(i,:), 'LineWidth',2.0, 'Color', 'magenta', 'LineStyle',':');
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

    plotConstr(tau*x_pos_lim, pos_lb, pos_ub, tau*xeq_pos, pos_eq, ax{1}, i);
    plotConstr(tau*x_vel_lim, vel_lb, vel_ub, tau*xeq_vel, vel_eq, ax{2}, i);
    plotConstr(tau*x_accel_lim, accel_lb, accel_ub, tau*xeq_accel, accel_eq, ax{3}, i);
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
