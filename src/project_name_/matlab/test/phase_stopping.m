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
K = 500;
D = 100;
gmp = GMP_nDoF(size(Pd_data,1), N_kernels, D, K, kernels_std_scaling);
offline_train_mse = gmp.train(train_method, Timed, Pd_data);

%% initialize DMP
a_z = D;
b_z = K/D;
can_clock_ptr = CanonicalClock();
% shape_attr_gat_ptr = SigmoidGatingFunction(1.0, 0.5);
% shape_attr_gat_ptr = LinGatingFunction(1.0, 0.05);
shape_attr_gat_ptr = LinGatingFunction(1.0, 0.01);
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

ad = 0.5;
[Time, y_data, x_data, dist_data] = simulateGMP_nDoF(gmp, P0, Pg, T, dt, false, ad);
[Time_rev, y_rev_data, x_rev_data, dist_rev_data] = simulateGMP_nDoF(gmp, P0, Pg, T, dt, true, ad);
[Time2, y_data2, x_data2, dist_data2] = simulateDMP(dmp, P0, Pg, T, dt, ad);

ticks_fs = 13;
labels_fs = 17;
legend_fs = 17;

%% Plot results
fig = figure();
fig.Position = [50 50 620 525];
ax = subplot(3,1,1); ax.FontSize = ticks_fs;
ax.Position = [0.1300 0.6283 0.7750 0.2727];
hold on;
plot(Timed, Pd_data, 'LineWidth',2.5 ,  'LineStyle','-', 'Color','green');
plot(Time, y_data, 'LineWidth',2.5 ,  'LineStyle','-', 'Color','blue');
plot(Time2, y_data2, 'LineWidth',2.5, 'LineStyle',':', 'Color','magenta');
plot(Time_rev, y_rev_data, 'LineWidth',2.5 ,  'LineStyle','--', 'Color',[0.85 0.33 0.1]);
legend({'training','novel', 'classical', 'reverse'}, 'interpreter','latex', 'fontsize',legend_fs, ...
    'Position',[0.0064 0.9345 0.9884 0.0576], 'Orientation','horizontal');
ylabel('position $y$', 'interpreter','latex', 'fontsize',labels_fs);
axis tight;
hold off;
ax = subplot(3,1,2); ax.FontSize = ticks_fs;
ax.Position = [0.1300 0.3430 0.7750 0.2157];
hold on;
plot(Time, x_data, 'LineWidth',2.5 ,  'LineStyle','-', 'Color','blue');
plot(Time2, x_data2, 'LineWidth',2.5 ,  'LineStyle',':', 'Color','magenta');
plot(Time_rev, x_rev_data, 'LineWidth',2.5 ,  'LineStyle','--', 'Color',[0.85 0.33 0.1]);
% legend({'$x$ - novel', '$x$ - classical', '$x$ - reverse'}, 'interpreter','latex', 'fontsize',15, 'Position',[0.6862 0.3627 0.2339 0.1004]);
ylabel('phase var. $x$', 'interpreter','latex', 'fontsize',labels_fs);
axis tight;
ax = subplot(3,1,3); ax.FontSize = ticks_fs;
ax.Position = [0.1300 0.1100 0.7750 0.1546];
plot(Time, dist_data, 'LineWidth',2.5 ,  'LineStyle','-', 'Color','red');
legend({'$d(t)$'}, 'interpreter','latex', 'fontsize',legend_fs);
xlabel('time [$s$]', 'interpreter','latex', 'fontsize',labels_fs);
ylabel('disturbance', 'interpreter','latex', 'fontsize',labels_fs);
axis tight;

%% ===================================================================

function d = disturbanceFun(t, T)

    a = [0.4 0.45 0.55 0.6];
    t1 = a(1)*T;
    t2 = a(2)*T;
    t3 = a(3)*T;
    t4 = a(4)*T;
    
    d_max = 10;
    
    if (t < t1), d = 0;
    elseif (t < t2), d = d_max*(t-t1)/(t2-t1);
    elseif (t < t3), d = d_max;
    elseif (t < t4), d = d_max + d_max*(t-t3)/(t3-t4);
    else, d = 0;
    end

end

function [Time, Y_data, x_data, dist_data] = simulateGMP_nDoF(gmp, y0, g, T, dt, reverse, ad)

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

d = 0;

iters = 0;
Time = [];
Y_data = [];
dY_data = [];
ddY_data = [];
x_data = [];
dist_data = [];

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
    dist_data = [dist_data d];

    %% DMP simulation
    y_c = 0; 
    z_c = 0.0;
    gmp.update(s, y, z, y_c, z_c);
    dy = gmp.getYdot();
    dz = gmp.getZdot();

    yc_dot = 0.0;
    ddy = gmp.getYddot(yc_dot);

    %% Update phase variable
    d = disturbanceFun(t,T);
    dx = dir * 1/tau * 1/(1 + ad*norm(d));

    %% Stopping criteria
    if (t>=1*t_end && norm(y-target)<5e-3 && norm(dy)<5e-3)
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




function [Time, Y_data, x_data, dist_data] = simulateDMP(dmp, y0, g, T, dt, ad)
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
d = 0;

t_end = T;
can_clock_ptr.setTau(t_end);

iters = 0;
Time = [];
Y_data = [];
dY_data = [];
ddY_data = [];
x_data = [];
dist_data = [];

for i=1:Dim, dmp{i}.setY0(y0(i)); end

%% simulate
while (true)

    %% data logging
    Time = [Time t];
    Y_data = [Y_data y];
    dY_data = [dY_data dy];  
    ddY_data = [ddY_data ddy];
    x_data = [x_data x];
    dist_data = [dist_data d];

    %% DMP simulation
    for i=1:Dim
        y_c = 0; 
        z_c = 0.0;
        dmp{i}.update(x, y(i), z(i), g(i), y_c, z_c);
        dy(i) = dmp{i}.getYdot();
        dz(i) = dmp{i}.getZdot();
    
        ddy(i) = dz(i)/dmp{i}.getTau();
    end

    d = disturbanceFun(t,T);
    
    %% Update phase variable
    can_clock_ptr.setTau( T*(1+ad*norm(d)) );
    dx = can_clock_ptr.getPhaseDot(x);
    
    

    %% Stopping criteria
    if (t>=1*t_end && norm(y-g)<5e-3 && norm(dy)<5e-3)
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


