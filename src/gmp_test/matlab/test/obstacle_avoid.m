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

%% initialize DMP
a_z = D;
b_z = K/D;
can_clock_ptr = CanonicalClock();
% shape_attr_gat_ptr = SigmoidGatingFunction(1.0, 0.5);
% shape_attr_gat_ptr = LinGatingFunction(1.0, 0.05);
shape_attr_gat_ptr = LinGatingFunction(1.0, 0.01);
dmp = cell(2,1);
for i=1:length(dmp)
    dmp{i} = DMP(N_kernels, a_z, b_z, can_clock_ptr, shape_attr_gat_ptr);
    offline_train_mse = dmp{i}.train(train_method, Timed, Pd_data(i,:), dPd_data(i,:), ddPd_data(i,:));
end


fig = figure;
fig.Position(3:4) = [495 417];
ax = axes('FontSize',14);
hold(ax, 'on');

plot(nan,nan, 'LineWidth',3.0 ,  'LineStyle','-', 'Color','blue', 'DisplayName','novel');
plot(nan,nan, 'LineWidth',3.0 ,  'LineStyle',':', 'Color','magenta', 'DisplayName','classical');
plot(nan, nan, 'LineWidth',3.0 ,  'LineStyle','-', 'Color',[0 0.7 0], 'DisplayName','training');
plot(nan, nan, 'LineWidth',3.0 ,  'LineStyle','--', 'Color',[0.85 0.33 0.1], 'DisplayName','reverse');
scatter(nan, nan, 'LineWidth',5.0 , 'LineWidth',2, 'Marker','o', ...
            'MarkerEdgeColor',[1 0 0], 'SizeData',200, 'DisplayName','obstacle');
        
legend({}, 'interpreter','latex', 'fontsize',18, 'Position',[0.6863 0.1697 0.2022 0.2603]);

plot(Pd_data(1,:), Pd_data(2,:), 'LineWidth',5.0 ,  'LineStyle','-', 'Color',[0 0.7 0], 'HandleVisibility','off');

ks = 1;
kt = 1; % spatial scale
% kt = 1.3; % temporal scale
P0 = Pd_data(:, 1);
Pgd = Pd_data(:, end);
Pg = P0 + ks*(Pgd - P0);
T = 4; %Timed(end) / kt;
dt = Ts;

obstacle = [0.366; 0.573];
gamma = 1000;
beta = 8;
[Time, P_data] = simulateGMP_nDoF(gmp, P0, Pg, T, dt, false, obstacle, gamma, beta);
[Time_rev, P_rev_data] = simulateGMP_nDoF(gmp, P0, Pg, T, dt, true, obstacle, gamma, beta);
[Time2, P_data2] = simulateDMP(dmp, P0, Pg, T, dt, obstacle, gamma, beta);

%% Plot results
plot(P_data(1,:), P_data(2,:), 'LineWidth',2.5 ,  'LineStyle','-', 'Color','blue', 'HandleVisibility','off');
plot(P_data2(1,:), P_data2(2,:), 'LineWidth',2.5, 'LineStyle',':', 'Color','magenta', 'HandleVisibility','off');
plot(P_rev_data(1,:), P_rev_data(2,:), 'LineWidth',2.5, 'LineStyle','--', 'Color',[0.85 0.33 0.1], 'HandleVisibility','off');
xlabel('x-axis [$m$]', 'interpreter','latex', 'fontsize',17);
ylabel('y-axis [$m$]', 'interpreter','latex', 'fontsize',17);

scatter(obstacle(1), obstacle(2), 'LineWidth',4, 'Marker','o', ...
    'MarkerEdgeColor',[1 0 0], 'SizeData',100, 'HandleVisibility','off');

axis tight;
hold off;


%% ===================================================================


function [Time, Y_data, dY_data, ddY_data] = simulateGMP_nDoF(gmp, y0, g, T, dt, reverse, obstacle, gamma, beta)

%% Simulates a dmp
% @param[in] gmp: Dim x 1 cell array, where each cell is a 1D GMP.
% @param[in] y0: Dim x 1 vector with the initial position..
% @param[in] g: Dim x 1 vector with the goal/target position.
% @param[in] T: Movement total time duration.
% @param[in] dt: Simulation timestep.
% @param[out] Time: 1 x N rowvector with simulation output timestamps.
% @param[out] Y_data: Dim x N matrix with simulation output positions.
% @param[out] dY_data: Dim x N matrix with simulation output velocities.
% @param[out] ddY_data: Dim x N matrix with simulation output accelerations.
%


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
    % x_data = [x_data x];

    %% DMP simulation
    y_c = 0; 
    z_c = 0;
    if (norm(dy) > 1e-3)
        phi = acos( dot((obstacle-y), dy) / (norm(obstacle-y) * norm(dy)) );
        z_c = gamma * [cos(phi); sin(phi)]*norm(dy)*phi*exp(-beta*phi);
    end
    gmp.update(s, y, z, y_c, z_c);
    dy = gmp.getYdot();
    dz = gmp.getZdot();

    yc_dot = 0.0;
    ddy = gmp.getYddot(yc_dot);

    %% Update phase variable
    dx = dir* 1/tau;

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




function [Time, Y_data, dY_data, ddY_data] = simulateDMP(dmp, y0, g, T, dt, obstacle, gamma, beta)
%% Simulates a dmp
% @param[in] dmp: Dim x 1 cell array, where each cell is a 1D DMP.
% @param[in] y0: Dim x 1 vector with the initial position..
% @param[in] g: Dim x 1 vector with the goal/target position.
% @param[in] T: Movement total time duration.
% @param[in] dt: Simulation timestep.
% @param[out] Time: 1 x N rowvector with simulation output timestamps.
% @param[out] Y_data: Dim x N matrix with simulation output positions.
% @param[out] dY_data: Dim x N matrix with simulation output velocities.
% @param[out] ddY_data: Dim x N matrix with simulation output accelerations.
%


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
    
    z_c = zeros(size(y));
    if (norm(dy) > 1e-3)
        phi = acos( dot((obstacle-y), dy) / (norm(obstacle-y) * norm(dy)) );
        z_c = gamma * [cos(phi); sin(phi)]*norm(dy)*phi*exp(-beta*phi);
    end

    %% DMP simulation
    for i=1:Dim
        y_c = 0;
        dmp{i}.update(x, y(i), z(i), g(i), y_c, z_c(i));
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


