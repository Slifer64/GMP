set_matlab_utils_path();

%% Load training data
% load('data/train_data.mat', 'Data');
load('data/fifthOrd_train_data.mat', 'Data');

Timed = Data.Time;
Pd_data = Data.Pos(1,:);
dPd_data = Data.Vel(1,:);
ddPd_data = Data.Accel(1,:);

Ts = Timed(2)-Timed(1);

train_method = 'LS';
N_kernels = 20;
%% initialize and train GMP
kernels_std_scaling = 1;
K = 500;
D = 100;
gmp = GMP(N_kernels, D, K, kernels_std_scaling);
% tic
offline_train_mse = gmp.train(train_method, Timed, Pd_data);
% offline_train_mse
% toc

%% initialize DMP
a_z = D;
b_z = K/D;
can_clock_ptr = CanonicalClock();
% shape_attr_gat_ptr = SigmoidGatingFunction(1.0, 0.5);
% shape_attr_gat_ptr = LinGatingFunction(1.0, 0.05);
shape_attr_gat_ptr = LinGatingFunction(1.0, 0.01);
dmp = DMP(N_kernels, a_z, b_z, can_clock_ptr, shape_attr_gat_ptr);
% tic
offline_train_mse = dmp.train(train_method, Timed, Pd_data, dPd_data, ddPd_data);
% offline_train_mse;
% toc

ks_data = [-2, -1.5, -0.5, 0.5 1.5 2];

fig = figure();
fig.Position(3:4) = [498 376];
ax = axes('Position',[0.0923 0.1464 0.8558 0.7519]);
ax.FontSize = 15;
hold(ax, 'on');

plot(nan,nan, 'LineWidth',3.0 ,  'LineStyle','-', 'Color','blue', 'DisplayName','novel');
plot(nan,nan, 'LineWidth',3.0 ,  'LineStyle',':', 'Color','magenta', 'DisplayName','classical');
plot(nan, nan, 'LineWidth',3.0 ,  'LineStyle','-.', 'Color',[0.85 0.33 0.1], 'DisplayName','reverse');
plot(nan, nan, 'LineWidth',3.0 ,  'LineStyle','-', 'Color',[0 0.7 0], 'DisplayName','training');
legend({}, 'interpreter','latex', 'fontsize',16, 'Orientation','horizontal', 'Position',[0.0309 0.9240 0.9462 0.0801]);

plot(Timed, Pd_data, 'LineWidth',3.0 ,  'LineStyle','-', 'Color',[0 0.7 0], 'HandleVisibility','off');

for k=1:length(ks_data)
    
    ks = ks_data(k);
    
    kt = 1; % spatial scale
    % kt = 1.3; % temporal scale
    P0 = Pd_data(1);
    Pgd = Pd_data(end);
    Pg = P0 + ks*(Pgd - P0);
    T = Timed(end) / kt;
    dt = Ts;
    
    [Time, P_data] = simulateGMP(gmp, P0, Pg, T, dt, false);
    [Time_rev, P_rev_data] = simulateGMP(gmp, P0, Pg, T, dt, true);
    [Time2, P_data2] = simulateDMP(dmp, P0, Pg, T, dt);

    %% Plot results
    plot(Time, P_data, 'LineWidth',3.0 ,  'LineStyle','-', 'Color','blue', 'HandleVisibility','off');
    plot(Time2, P_data2, 'LineWidth',3.0, 'LineStyle',':', 'Color','magenta', 'HandleVisibility','off');
    plot(Time_rev, P_rev_data, 'LineWidth',3.0, 'LineStyle','-.', 'Color',[0.85 0.33 0.1], 'HandleVisibility','off');
    xlabel('time [$s$]', 'interpreter','latex', 'fontsize',17, 'Position',[1.0025 -2.3548 -1]);

end

% plot(ax.XLim, [Pgd Pgd], 'LineWidth',3.0, 'LineStyle','-', 'Color',[0 0 1 0.3], 'HandleVisibility','off');

axis tight;
hold off;


%% ================================================================

function [Time, Y_data, dY_data, ddY_data] = simulateGMP(gmp, y0, g, T, dt, reverse)
%% Simulates a gmp


%% set initial values
if (~iscell(gmp)), gmp = {gmp}; end
Dim = length(gmp);
y = y0;
dy = zeros(Dim,1);
ddy = zeros(Dim,1);
z = zeros(Dim,1);
dz = zeros(Dim,1);
t = 0.0;

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

for i=1:Dim
    gmp{i}.setY0(y0(i));
    gmp{i}.setGoal(g(i));
end

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
    % x_data = [x_data x];

    %% DMP simulation
    for i=1:Dim
        y_c = 0.0;
        z_c = 0.0;
        gmp{i}.update(s, y(i), z(i), y_c, z_c);
        dy(i) = gmp{i}.getYdot();
        dz(i) = gmp{i}.getZdot();
        
        yc_dot = 0.0;
        ddy(i) = gmp{i}.getYddot(yc_dot);
    end

    %% Update phase variable
    x_dot = dir*1/tau;

    %% Stopping criteria
    if (t>=t_end) % && norm(y-g)<5e-3 && norm(dy)<5e-3)
        break;
    end

    %% Numerical integration
    iters = iters + 1;
    t = t + dt;
    x = x + x_dot*dt;
    y = y + dy*dt;
    z = z + dz*dt;

    s = [x; x_dot];
    
end


end



function [Time, Y_data, dY_data, ddY_data] = simulateDMP(dmp, y0, g, T, dt)
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
    % x_data = [x_data x];

    %% DMP simulation
    for i=1:Dim
        y_c = 0.0;
        z_c = 0.0;
        dmp{i}.update(x, y(i), z(i), g(i), y_c, z_c);
        dy(i) = dmp{i}.getYdot();
        dz(i) = dmp{i}.getZdot();
    
        ddy(i) = dz(i)/dmp{i}.getTau();
    end

    %% Update phase variable
    dx = can_clock_ptr.getPhaseDot(x);

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

end


end



