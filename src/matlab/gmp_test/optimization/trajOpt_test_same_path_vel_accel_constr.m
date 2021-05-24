clc;
close all;
clear;

rng(0);

set_matlab_utils_path();

%% Load training data
load('data/traj_opt_train_data.mat', 'Data');

ind = [1 2 3];

Timed = Data.Time;
Pd_data = Data.Pos(ind,:);
dPd_data = Data.Vel(ind,:);
ddPd_data = Data.Accel(ind,:);

Ts = Timed(2)-Timed(1);

%% initialize and train GMP
train_method = 'LS';
N_kernels = 30;
kernels_std_scaling = 1;
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

kt = 0.5; % temporal scaling
ks = diag([2 2 2]); % spatial scaling
tau = taud/kt;
y0 = yd0;
yg = ks*(ygd - yd0) + y0;
ydot_0 = 0;
gmp.setY0(y0);
gmp.setGoal(yg);

Zero = zeros(n_dof,1);

% =======  Position constraints  =======
pos_up_lim = [3; 3; 3];
pos_low_lim = -pos_up_lim;
% pos_constr = [GMPConstr(0,y0,'='); GMPConstr(1,yg,'=')];
xi = 0:0.02:1;
pos_eq_constr = [];%repmat(GMPConstr(0,0,'='), length(xi), 1);
% for j=1:length(xi)
%     pos_eq_constr(j) = GMPConstr(xi(j),gmp.getYd(xi(j)),'=');
% end
pos_constr = [pos_eq_constr; GMP_Opt.lowerBoundConstr(xi, pos_low_lim(ind)); GMP_Opt.upperBoundConstr(xi, pos_up_lim(ind))];
% =======  Velocity constraints  =======
vel_constr = [];%[GMPConstr(0,Zero,'='); GMPConstr(1,Zero,'=')];
xi = 0:0.02:1;
vel_constr = [vel_constr; GMP_Opt.thresholdConstr(xi, 0.3*ones(n_dof,1))];
% =======  Acceleration constraints  =======
accel_constr = [GMPConstr(0,Zero,'='); GMPConstr(1,Zero,'=')];
xi = 0:0.02:1;
accel_constr = [accel_constr; GMP_Opt.thresholdConstr(xi, 0.6*ones(n_dof,1))];

tic
gmp_opt = GMP_Opt(gmp);
gmp_opt.setOptions(true, false, false, 1, 1, 1);
% gmp_opt.constrOpt(200, tau, pos_constr, vel_constr, accel_constr);
gmp_opt.constrOpt(200, tau, pos_constr, vel_constr, []);
fprintf([gmp_opt.getExitMsg() '\n']);

toc

T_sc = gmp.getScaling();

%% calculate scaled demo trajectory
Time2 = Timed / kt;
P_data2 = T_sc*(Pd_data - Pd_data(:,1)) + y0;
dP_data2 = kt*T_sc*dPd_data;
ddP_data2 = (kt^2)*T_sc*ddPd_data;


%% simulation

Time = [];
P_data = [];
dP_data = [];
ddP_data = [];

p = y0;
p_dot = ydot_0;
p_ddot = 0;

t = 0;
dt = 0.01;

while (true)

    x = t/tau;
    x_dot = 1/tau;
    p = gmp.getYd(x);
    p_dot = gmp.getYdDot(x, x_dot);
    p_ddot = gmp.getYdDDot(x, x_dot, 0);

    Time = [Time t];
    P_data = [P_data p];
    dP_data = [dP_data p_dot];
    ddP_data = [ddP_data p_ddot];

    t = t + dt;

    if (t >= 1.05*tau), break; end

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

    plotConstr(pos_constr, ax{1}, tau, i);
    plotConstr(vel_constr, ax{2}, tau, i);
    plotConstr(accel_constr, ax{3}, tau, i);
end

%% =======================================================
%% =======================================================

function plotConstr(constr, ax, tau, i_dof)

    hold(ax, 'on');

    eq_sc = scatter(nan, nan, 'MarkerEdgeColor',[0 0.7 0], 'Marker','*', ...
            'MarkerEdgeAlpha',0.6, 'LineWidth',2, 'SizeData', 100, 'Parent',ax);

    less_sc = scatter(nan, nan, 'MarkerEdgeColor',[1 0 0], 'Marker','v', ...
            'MarkerEdgeAlpha',0.3, 'LineWidth',1, 'SizeData', 100, 'Parent',ax);

    gr_sc = scatter(nan, nan, 'MarkerEdgeColor',[1 0 0], 'Marker','^', ...
            'MarkerEdgeAlpha',0.3, 'LineWidth',1, 'SizeData', 100, 'Parent',ax);

    for i=1:length(constr)
        t = constr(i).x * tau;
        value = constr(i).value(i_dof);
        const_type = constr(i).type;

        if (const_type == '=')
           eq_sc.XData = [eq_sc.XData t];
           eq_sc.YData = [eq_sc.YData value];
        elseif (const_type == '<')
            less_sc.XData = [less_sc.XData t];
            less_sc.YData = [less_sc.YData value];
        elseif (const_type == '>')
            gr_sc.XData = [gr_sc.XData t];
            gr_sc.YData = [gr_sc.YData value];
        end

    end

end
