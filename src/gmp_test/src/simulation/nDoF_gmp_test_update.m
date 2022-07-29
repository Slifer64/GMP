%% Simulates a GMP for encoding and executing a Cartesian position trajectory.
%  Works similarly for an n-DoF trajectory as well.

clc;
close all;
clear;

%% =============  includes...  =============
import_gmp_lib();
import_io_lib();


%% =============  Set GMP params  =============
train_method = 'LS'; % {'LS', 'LWR'}
N_kernels = 25;  % number or kernels (Gaussians)
kernels_std_scaling = 1.5; % scaling of the width of each Gaussian. The greater 
                           % the scaling the more the overlapping). A vaule 
                           % in [1 1.5] is usually good.

train_filename = 'data/pos_data.bin'; % file containing the training data

%% Choose scaling type for the generalization
scale_type = 'none'; % {'prop', 'rot_min', 'rot_wb', 'none'}, use 'prop' if n-DoFs > 3
wb_normal = [0; 0; 1]; % work-bench normal (for use with "rot_wb")

% %% Set scaling for the executed trajectory (for testing purposes...)
% spat_s = [2; 3; 1.1]; % spatial scaling
% temp_s = 1.2; % temporal scaling: execute the trajectory 1.2 times faster than the demo

%% For loading/writing the gmp from/to a file.
read_gmp_from_file = false;
write_gmp_to_file = false;

gmp_filename = 'data/gmp_pos.bin'; % location of the DMP model

%% =============  Load train data  =============
fid = FileIO(train_filename, FileIO.in);
Timed = fid.read('Timed');
Pd_data = fid.read('Pd_data');
dPd_data = fid.read('dPd_data');
ddPd_data = fid.read('ddPd_data');
% % Or you can just do
% data = fid.readAll();
% Timed = data.Timed;
% Pd_data = data.Pd_data;
% dPd_data = data.dPd_data;
% ddPd_data = data.ddPd_data;
fid.close();

Ts = Timed(2) - Timed(1);

%% =============  Create/Train GMP  =============

n_dof = size(Pd_data, 1); % number of DoFs

if (read_gmp_from_file)

    gmp = GMP();
    gmp_.read(gmp, gmp_filename, '');
    disp('Loaded GMP from file!');

else

    %% initialize and train GMP
    gmp = GMP(n_dof, N_kernels, kernels_std_scaling);
    t_start = tic;
    offline_train_mse = gmp.train(train_method, Timed/Timed(end), Pd_data);
    offline_train_mse
    toc(t_start)

end

% % Optionally, readjust the number of kernels and their widths
% gmp.autoRetrain(50, 2, 200, 'LWR');

% set scaling method for generalizing the DMP to new targets/initial positions
% You can also add your own scaling method, by implementing a similar class
% that inherits from @TrajScale.
if ( strcmpi(scale_type, 'prop') )
    traj_sc = TrajScale_Prop(n_dof);
elseif ( strcmpi(scale_type, 'rot_min') )
    traj_sc = TrajScale_Rot_min();
elseif ( strcmpi(scale_type, 'rot_wb') )
    traj_sc = TrajScale_Rot_wb();
    traj_sc.setWorkBenchNormal(wb_normal); % set also the workbench normal
elseif ( strcmpi(scale_type, 'none') )
    traj_sc = TrajScale_None(n_dof);
else
    error(['Unsupported scale type ''' scale_type '''...\n']);
end

gmp.setScaleMethod(traj_sc);


if (write_gmp_to_file)
    gmp_.write(gmp, gmp_filename, ''); 
    disp('Wrote GMP from file!');
end

%% DMP simulation
disp('GMP simulation...');
t_start = tic;

%% Initial/Final values
P0d = Pd_data(:,1);   % Initial demo position
Pgd = Pd_data(:,end); % Target demo position
P0 = P0d; % set initial position for execution (for simplicity lets leave it the same as the demo)
% Pg = [1.2; 1.8; 1.36];  % set target position for execution
Pg = [0.8; 0.9; -0.2];  % set target position for execution
% Pg = spat_s.*(Pgd - P0) + P0;
T = 8.33; % set the time duration of the executed motion
% T = Timed(end) / temp_s;


dt = Ts; % time step for numerical integration

%% Execute the DMP
[Time, P_data, dP_data, ddP_data] = simulateGMP(gmp, P0, Pg, T, dt);
toc(t_start)

%% This is the groundtruth trajectory that should be produced
Ks = gmp.getScaling(); % spatial scaling 
temp_s = Timed(end) / T; % temporal scaling 
Timed2 = Timed / temp_s;
Pd2_data = Ks*( Pd_data-P0d ) + P0;
dPd2_data = Ks*dPd_data*temp_s;
ddPd2_data = Ks*ddPd_data*temp_s^2;

%% Plot results
for i=1:3
    figure;
    subplot(3,1,1);
    hold on;
    plot(Time(end), Pg(i), 'LineWidth',4, 'Marker','x', 'MarkerSize',10, 'LineStyle','none', 'Color','red', 'HandleVisibility','off');
    plot(Time, P_data(i,:), 'LineWidth',2.0 , 'Color','blue', 'DisplayName', 'DMP');
    plot(Timed2, Pd2_data(i,:), 'LineWidth',2.0, 'LineStyle',':', 'Color','magenta', 'DisplayName','ground-truth');
    plot(Timed*Time(end)/Timed(end), Pd_data(i,:), 'LineWidth',2.0, 'LineStyle','-.' , 'Color','green', 'DisplayName', 'demo');
    ylabel('pos [$m$]', 'interpreter','latex', 'fontsize',15);
    title(['DoF $' num2str(i) '$'], 'interpreter','latex', 'fontsize',18);
    legend({}, 'interpreter','latex', 'fontsize',15);

    axis tight;
    hold off;

    subplot(3,1,2);
    hold on;
    plot(Time, dP_data(i,:), 'LineWidth',2.0, 'Color','blue');
    plot(Timed2, dPd2_data(i,:), 'LineWidth',2.0, 'LineStyle',':', 'Color','magenta');
    ylabel('vel [$m/s$]', 'interpreter','latex', 'fontsize',15);
    axis tight;
    hold off;

    subplot(3,1,3);
    hold on;
    plot(Time, ddP_data(i,:), 'LineWidth',2.0, 'Color','blue');
    plot(Timed2, ddPd2_data(i,:), 'LineWidth',2.0, 'LineStyle',':', 'Color','magenta');
    ylabel('accel [$m/s^2$]', 'interpreter','latex', 'fontsize',15);
    axis tight;
    hold off;
end


fig = figure;
fig.Position(3:4) = [686 616];
hold on;
plot3(P_data(1,:), P_data(2,:), P_data(3,:), 'LineWidth',2, 'LineStyle','-', 'Color','magenta', 'DisplayName', 'DMP');
plot3(Pd2_data(1,:), Pd2_data(2,:), Pd2_data(3,:), 'LineWidth',2, 'LineStyle',':', 'Color','blue', 'DisplayName','ground-truth');
plot3(Pd_data(1,:), Pd_data(2,:), Pd_data(3,:), 'LineWidth',2, 'LineStyle',':', 'Color',[0 0.7 0], 'DisplayName','demo');
plot3(P0(1), P0(2), P0(3), 'LineWidth',4, 'Marker','o', 'MarkerSize',10, 'LineStyle','none', 'Color','green', 'DisplayName', '$p_0$');
plot3(Pg(1), Pg(2), Pg(3), 'LineWidth',4, 'Marker','x', 'MarkerSize',10, 'LineStyle','none', 'Color','red', 'DisplayName', 'target');
plot3(Pgd(1), Pgd(2), Pgd(3), 'LineWidth',4, 'Marker','x', 'MarkerSize',10, 'LineStyle','none', 'Color','magenta', 'DisplayName', 'demo target');
legend({}, 'interpreter','latex', 'fontsize',15, 'Position',[0.1887 0.7924 0.2204 0.1711]);
xlabel('$x$', 'interpreter','latex', 'fontsize',15);
ylabel('$y$', 'interpreter','latex', 'fontsize',15);
zlabel('$z$', 'interpreter','latex', 'fontsize',15);
view(-7.8, 30.92);
axis tight;
grid on;
hold off;


%% ============================================================
%% ============================================================


function [Time, Y_data, dY_data, ddY_data] = simulateGMP(gmp, y0, g, T, dt)
%% Simulates a dmp
% @param[in] gmp: object of type @GMP.
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
Dim = gmp.numOfDoFs(); % number of DoFs

y = y0; % position
dy = zeros(Dim,1); % velocity
ddy = zeros(Dim,1); % acceleration
O_ndof = zeros(Dim, 1);

t = 0.0;

% dz = zeros(Dim,1);
% z = zeros(Dim,1);

t_end = T;
tau = t_end;

% phase variable, from 0 to 1
x = 0.0;
x_dot = 1/tau;
x_ddot = 0; % since x_dot is constant here
% s = GMP_phase(x, x_dot, 0);

iters = 0;

% data to log
Time = [];
Y_data = [];
dY_data = [];
ddY_data = [];
x_data = [];

% set initial and final position
if (gmp.traj_sc.getScaleType() ~= TrajScale.NONE)
    gmp.setY0(y0);
    gmp.setGoal(g);
else
    W0 = gmp.W;
    gmp_up = GMP_Update(gmp);
    gmp_up.syncUpdate(false);
%     gmp_up.initExpSigmaw(0.01);
%     gmp_up.initSigmaWfromMsr(0:0.01:1);
%     gmp_up.initSigmaWfromVelMsr(0:0.01:1, tau);
    gmp_up.initSigmaWfromAccelMsr(0:0.01:1, tau);
    gmp_up.plotWeightsCovariance(); %pause;
    gmp_up.updatePos(0, y0, 1e-6);
    gmp_up.updateVel(0, x_dot, O_ndof, 1e-4);
    gmp_up.updateAccel(0, x_dot, x_ddot, O_ndof, 1e-4);
    gmp_up.updatePos(1, g, 1e-6);
    gmp_up.updateVel(1, x_dot, 2*[0.05; 0.05; 0.05], 1e-6);
    gmp_up.updateAccel(1, x_dot, x_ddot, O_ndof, 1e-4);
    gmp_up.updateNow();
    %gmp.updateGoal(g, 1/tau, GMP_phase(x, x_dot, x_ddot), y, dy);
end

count = 0;
g0 = g;
g_data = [];
W2_data = [];


%% simulate
while (true)

    %% data logging
    Time = [Time t];
    Y_data = [Y_data y];
    dY_data = [dY_data dy];
    ddY_data = [ddY_data ddy];
    % x_data = [x_data x];
    
    g_data = [g_data g];
    W2_data = [W2_data gmp.W(2,:)'];

    %% DMP simulation
    y_x = gmp.getYd(x);
    dy_x = gmp.getYdDot(x, x_dot);
    ddy_x = gmp.getYdDDot(x, x_dot, x_ddot);
    
    K = 300; % set the DMP stiffness
    D = 60; % set the DMP damping
    
    external_signal = 0; % optionally add some external signal
    
    % Track it using a 2nd order dynamical system. This is actually the DMP. 
    ddy = ddy_x + D*(dy_x - dy) + K*(y_x - y) + external_signal;

    % if external_signal == 0 you can obviously set directly:
%     y = y_x;
%     dy = dy_x;
%     ddy = ddy_x;
    
    %% Stopping criteria
    if (t>=1.0*t_end) % && norm(y-g)<1e-3 && norm(dy)<5e-3)
        break;
    end
    
    count = count + 1;
    
    %% Numerical integration
    iters = iters + 1;
    t = t + dt;
    x = x + x_dot*dt;
    x_dot = x_dot + x_ddot*dt;
    y = y + dy*dt;
    dy = dy + ddy*dt;
%     z = z + dz*dt;

%     s.x = x;

end

% figure;
% plot(Time, g_data)
% 
% figure;
% plot(W2_data')

end
