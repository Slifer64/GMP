set_matlab_utils_path();

%% Load training data
% load('data/train_data.mat', 'Data');
load('data/train_data1.mat', 'Data');

ind = [1 3];

Timed = Data.Time;
Pd_data = Data.Pos(ind,:);
dPd_data = Data.Vel(ind,:);
ddPd_data = Data.Accel(ind,:);

Ts = Timed(2)-Timed(1);

train_method = 'LS';
N_kernels = 30;
%% initialize and train GMP
kernels_std_scaling = 1;
K = 500;
D = 100;
gmp = GMP_nDoF(size(Pd_data,1), N_kernels, D, K, kernels_std_scaling);
% tic
offline_train_mse = gmp.train(train_method, Timed, Pd_data);
% offline_train_mse
% toc

ks = 1;
kt = 1; % spatial scale
% kt = 1.3; % temporal scale
P0 = Pd_data(:, 1);
Pgd = Pd_data(:, end);
Pg = P0 + ks*(Pgd - P0);
T = Timed(end) / kt;
dt = Ts;

[Time, P_data] = simulateGMP_nDoF(gmp, P0, Pg, T, dt, false);
[Time2, P_data2] = simulateGMP_nDoF(gmp, P0, Pg, T, dt, true);

%% Plot results
figure;
ax = axes();
hold on
plot(P_data(1,:), P_data(2,:), 'LineWidth',3 ,  'LineStyle','-', 'Color','cyan');
plot(P_data2(1,:), P_data2(2,:), 'LineWidth',4, 'LineStyle',':', 'Color','magenta');
scatter(P0(1), P0(2), 'LineWidth',4, 'Marker','o', 'MarkerEdgeColor',[0 1 0], 'SizeData',200);
scatter(Pg(1), Pg(2), 'LineWidth',4, 'Marker','x', 'MarkerEdgeColor',[1 0 0], 'SizeData',200);
xlabel('x-axis [$m$]', 'interpreter','latex', 'fontsize',15);
ylabel('y-axis [$m$]', 'interpreter','latex', 'fontsize',15);
legend({'forward','reverse','initial','target'}, 'interpreter','latex', 'fontsize',15, 'Position',[0.6840 0.2293 0.1891 0.2256]);
axis tight;
ax.XLim = ax.XLim + [-0.05 0.05];
ax.YLim = ax.YLim + [-0.05 0.05];
hold off;


%% ===================================================================

function [Time, Y_data, dY_data, ddY_data] = simulateGMP_nDoF(gmp, y0, g, T, dt, reverse)
%% Simulates a gmp

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

dir = 1;
if (reverse)
    dir = -1;
    y = g;
    x = 1;
end

%% simulate
while (true)

    %% data logging
    Time = [Time t];
    Y_data = [Y_data y];
    dY_data = [dY_data dy];  
    ddY_data = [ddY_data ddy];
    % x_data = [x_data x];

    %% DMP simulation
    y_c = 0.0;
    z_c = 0.0;
    gmp.update(s, y, z, y_c, z_c);
    dy = gmp.getYdot();
    dz = gmp.getZdot();

    yc_dot = 0.0;
    ddy = gmp.getYddot(yc_dot);

    %% Update phase variable
    dx = 1/tau * dir;

    %% Stopping criteria
    if (t>=1*t_end) % && norm(y-g)<5e-3 && norm(dy)<5e-3)
        break;
    end

    %% Numerical integration
    iters = iters + 1;
    t = t + dt;
    x = x + dx*dt;
    y = y + dy*dt;
    z = z + dz*dt;

    s = [x; dx];
    
end


end




