clear;
clc;
% close all;


set_matlab_utils_path();

%% Load training data
load('data/fifthOrd_train_data.mat', 'Data');

Timed = Data.Time;
Pd_data = Data.Pos;
dPd_data = Data.Vel;
ddPd_data = Data.Accel;

Ts = Timed(2)-Timed(1);

train_method = 'LS';
N_kernels = 30;
%% initialize and train GMP
kernels_std_scaling = 1;
K = 200;
D = 30;
gmp = GMP_nDoF(size(Pd_data,1), N_kernels, D, K, kernels_std_scaling);
offline_train_mse = gmp.train(train_method, Timed, Pd_data);

%% initialize DMP
a_z = D;
b_z = K/D;
can_clock_ptr = CanonicalClock();
% shape_attr_gat_ptr = SigmoidGatingFunction(1.0, 0.5);
% shape_attr_gat_ptr = LinGatingFunction(1.0, 0.01);
% shape_attr_gat_ptr = ExpGatingFunction(1.0, 0.01);
shape_attr_gat_ptr = LinGatingFunction(1.0, 1.0);
dmp = DMP(N_kernels, a_z, b_z, can_clock_ptr, shape_attr_gat_ptr);
offline_train_mse = dmp.train(train_method, Timed, Pd_data(1,:), dPd_data(1,:), ddPd_data(1,:));

ks = 1;
kt = 1; % spatial scale
% kt = 1.3; % temporal scale
P0 = Pd_data(:, 1);
Pgd = Pd_data(:, end);
Pg = P0 + ks*(Pgd - P0);
T = Timed(end) / kt;
dt = 0.001;

new_goal = struct('t',0.6*T, 'g',Pg+0.4);
[Time, y_data, x_data] = simulateGMP_nDoF(gmp, P0, Pg, T, dt, false, new_goal);
[Time_rev0, y_rev0_data, x_rev0_data] = simulateGMP_nDoF(gmp, P0, new_goal.g, T, dt, true, struct('t',-1, 'g',[]));
[Time_rev, y_rev_data, x_rev_data] = simulateGMP_nDoF(gmp, P0, Pg, T, dt, true, struct('t',0.7*T, 'g',P0-0.5));
[Time2, y_data2, x_data2] = simulateDMP(dmp, P0, Pg, T, dt, new_goal);

ticks_fs = 13;
labels_fs = 17;
legend_fs = 17;

%% Plot results
figure('Position',[100 100 620 525]);
% --------------------------------------------
ax = subplot(2,1,1);
ax.set('FontSize',ticks_fs);
hold on;
% plot(Timed, Pd_data, 'LineWidth',2.5 ,  'LineStyle','-', 'Color','green', 'DisplayName','training');
plot(Time, y_data, 'LineWidth',2.5 ,  'LineStyle','-', 'Color','blue', 'DisplayName', 'novel');
plot(Time_rev0, y_rev0_data, 'LineWidth',2.5 ,  'LineStyle',':', 'Color',[0.85 0.33 0.1], 'DisplayName', 'novel:reverse');
plot(Time2, y_data2, 'LineWidth',2.5, 'LineStyle',':', 'Color','magenta', 'DisplayName', 'classical');
scatter(Timed(end), Pd_data(end), 'LineWidth',3.0, 'Marker','o', ...
            'MarkerEdgeColor',[0 1 1], 'SizeData',150, 'DisplayName','initial $g$');
scatter(Time(end), y_data(end), 'LineWidth',3.0, 'Marker','o', ...
    'MarkerEdgeColor',[1 0 0], 'SizeData',150, 'DisplayName','new $g$');
legend({}, 'interpreter','latex', 'fontsize',legend_fs, 'Orientation','horizontal', 'Position',[0.0047 0.9197 0.9872 0.0604]);
ylabel('position $y$', 'interpreter','latex', 'fontsize',labels_fs);
% title({['novel-DMP: Stiffness = ' num2str(K) ', Damping = ' num2str(D)];
%     ['DMP: Stiffness = ' num2str(K/T^2) ', Damping = ' num2str(D/T)]}, 'fontsize',16);
ax.Position = [0.1300 0.5476 0.7750 0.3412];
axis tight;
% --------------------------------------------
ax = subplot(2,1,2);
ax.set('FontSize',ticks_fs, 'NextPlot','add');
plot(Timed, fliplr(Pd_data), 'LineWidth',2.5 ,  'LineStyle','-', 'Color',[0 0.8 0], 'DisplayName','training - reverse');
plot(Time_rev, y_rev_data, 'LineWidth',2.5 ,  'LineStyle','--', 'Color',[0.85 0.33 0.1], 'DisplayName','novel - reverse');
scatter(Timed(end), Pd_data(1), 'LineWidth',3.0, 'Marker','o', ...
            'MarkerEdgeColor',[0 1 1], 'SizeData',150, 'DisplayName','initial $y_0$');
scatter(Time_rev(end), y_rev_data(end), 'LineWidth',3.0, 'Marker','o', ...
    'MarkerEdgeColor',[1 0 0], 'SizeData',150, 'DisplayName','new $y_0$');
legend({}, 'interpreter','latex', 'fontsize',legend_fs, 'Position',[0.1381 0.1188 0.3388 0.2177]);
ylabel('position $y$', 'interpreter','latex', 'fontsize',labels_fs);
xlabel('time [$s$]', 'interpreter','latex', 'fontsize',labels_fs);
axis tight;

%% ===================================================================

function [Time, Y_data, x_data] = simulateGMP_nDoF(gmp, y0, g, T, dt, reverse, new_goal)

%% set initial values
Dim = gmp.length();
ddy = zeros(Dim,1);
dy = zeros(Dim,1);
y = y0;
t = 0.0;
dz = zeros(Dim,1);
z = zeros(Dim,1);

t_end = T;
tau = t_end;
x = 0.0;
x_dot = 1/tau;
s = [x; x_dot];

iters = 0;
Time = [];
Y_data = [];
dY_data = [];
ddY_data = [];
x_data = [];

gmp.setY0(y0);
gmp.setGoal(g);

target = g;
dir = 1;
if (reverse)
    dir = -1;
    y = g;
    x = 1;
    target = y0;
end

%% simulate
while (true)

    %% data logging
    Time = [Time t];
    Y_data = [Y_data y];
    dY_data = [dY_data dy];  
    ddY_data = [ddY_data ddy];
    x_data = [x_data x];
    
    if (abs(t-new_goal.t) <= dt/2)
        if (dir>0), gmp.setGoal(new_goal.g);
        else, gmp.setY0(new_goal.g); 
        end
        target = new_goal.g;
    end

    %% DMP simulation
    y_c = 0; 
    z_c = 0.0;
    gmp.update(s, y, z, y_c, z_c);
    dy = gmp.getYdot();
    dz = gmp.getZdot();

    yc_dot = 0.0;
    ddy = gmp.getYddot(yc_dot);

    %% Update phase variable
    dx = dir * 1/tau;

    %% Stopping criteria
    if (t>=1*t_end && norm(y-target)<5e-2 && norm(dy)<5e-2)
        break;
    end

%     if (t >= t_end)
%         dir
%         pos_err = norm(y-target)
%         vel = norm(dy)
%         pause
%     end

    %% Numerical integration
    iters = iters + 1;
    t = t + dt;
    x = x + dx*dt;
    y = y + dy*dt;
    z = z + dz*dt;

    if (x>=1)
        x = 1;
        dx = 0;
    end
    
    if (x<=0)
        x = 0;
        dx = 0;
    end
    
    s = [x; dx];
    
    
%     if (t >=t_end && dir < 0)
%         
%         yd = gmp.getYd(x);
%         [yd y]
%         
%         yd_dot = gmp.getYdDot(x, dx);
%         [yd_dot dy]
%         pause
%         
%     end
    
end


end



function [Time, Y_data, x_data] = simulateDMP(dmp, y0, g, T, dt, new_goal)
%% Simulates a dmp


%% set initial values
if (~iscell(dmp)), dmp = {dmp}; end
can_clock_ptr = dmp{1}.can_clock_ptr;
Dim = length(dmp);
x = 0.0;
dx = 0.0;
ddy = zeros(Dim,1);
dy = zeros(Dim,1);
y = y0;
t = 0.0;
dz = zeros(Dim,1);
z = zeros(Dim,1);

t_end = T;
can_clock_ptr.setTau(t_end);

iters = 0;
Time = [];
Y_data = [];
dY_data = [];
ddY_data = [];
x_data = [];

for i=1:Dim, dmp{i}.setY0(y0(i)); end

%% simulate
while (true)

    %% data logging
    Time = [Time t];
    Y_data = [Y_data y];
    dY_data = [dY_data dy];  
    ddY_data = [ddY_data ddy];
    x_data = [x_data x];
    
    if (abs(t-new_goal.t) <= dt/2)
        g = new_goal.g;
    end

    %% DMP simulation
    for i=1:Dim
        y_c = 0; 
        z_c = 0.0;
        dmp{i}.update(x, y(i), z(i), g(i), y_c, z_c);
        dy(i) = dmp{i}.getYdot();
        dz(i) = dmp{i}.getZdot();
    
        ddy(i) = dz(i)/dmp{i}.getTau();
    end
    
    %% Update phase variable
    dx = can_clock_ptr.getPhaseDot(x);
    
    %% Stopping criteria
    if (t>=1*t_end && norm(y-g)<5e-2 && norm(dy)<5e-2)
        break;
    end

    %% Numerical integration
    iters = iters + 1;
    t = t + dt;
    x = x + dx*dt;
    y = y + dy*dt;
    z = z + dz*dt;

end


end


