clc;
close all;
clear;

%% =============  includes...  =============
import_gmp_lib();
import_io_lib();

%% =============  Load params  =============
train_method = 'LS';  % {'LS', 'LWR'}
N_kernels = 25; % number or kernels (Gaussians)
kernels_std_scaling = 1.5; % scaling of the width of each Gaussian. The greater 
                           % the scaling the more the overlapping). A vaule 
                           % in [1 1.5] is usually good.

train_filename = 'data/orient_data.bin'; % file containing the training data

%% Choose scaling type for the generalization
scale_type = 'rot_wb'; % {"prop", "rot_min", "rot_wb"}
wb_normal = [0; 0; 1]; % work-bench normal (for use with "rot_wb")

% %% Set scaling for the executed trajectory (for testing purposes...)
% spat_s = [1.2 1.3 1.1]; % spatial scaling
% temp_s = 1.2; % temporal scaling: execute the trajectory 1.2 times faster than the demo

read_gmp_from_file = false;
write_gmp_to_file = false;

gmp_filename = 'data/gmp_orient.bin'; % location of the DMP model

%% Load training data
fid = FileIO(train_filename, FileIO.in );
Timed = fid.read('Timed');
Qd_data = fid.read('Qd_data');
vRotd_data = fid.read('vRotd_data');
dvRotd_data = fid.read('dvRotd_data');
% % Or you can just do
% data = fid.readAll();
% ...
fid.close();

Ts = Timed(2)-Timed(1);

%% Pointer to function for simulating the DMP
% Simulate the DMP in the Cartesian orientation space:
simulateGMPo = @simulateGMPo_in_Cart_space;

% Or you can simulate it in the quaternion-logarithm space:
% simulateGMPo = @simulateGMPo_in_log_space;

%% Load / train DMP
if (read_gmp_from_file)

    gmp_o = GMPo();
    gmp_.read(gmp_o, gmp_filename, '');
    disp('Loaded GMP from file!');

else

    %% initialize and train GMP
    gmp_o = GMPo(N_kernels, kernels_std_scaling);
    tic
    offline_train_mse = gmp_o.train(train_method, Timed/Timed(end), Qd_data);
    offline_train_mse
    toc
    
    % set scaling method for generalizing the DMP to new targets/initial
    % orientations
    if ( strcmpi(scale_type, 'prop') )
       traj_sc = TrajScale_Prop(3);
    elseif ( strcmpi(scale_type, 'rot_min') )
        traj_sc = TrajScale_Rot_min();
    elseif ( strcmpi(scale_type, 'rot_wb') )
        traj_sc = TrajScale_Rot_wb();
        traj_sc.setWorkBenchNormal(wb_normal); % set also the workbench normal
    else
        error(['Unsupported scale method ''' scale_type '''...']);
    end

    gmp_o.setScaleMethod(traj_sc);

end

if (write_gmp_to_file)
    gmp_.write(gmp_o, 'gmp_o.bin', '');
    disp('Wrote GMP from file!');
end

%% DMP simulation
disp('GMP simulation...');
t_start = tic;

%% Initial/Final values
Qd0 = Qd_data(:,1);  % Initial demo orientation
Qgd = Qd_data(:,end); % Target demo orientation
Q0 = Qd0; % set initial position for execution (for simplicity lets leave it the same as the demo)
Qg = [0.7361; 0.3551; 0.4142; -0.4007]; % set target orientation for execution
% e0 = diag(spat_s)*gmp_.quatLogDiff(Qgd,Qd0);
% Qg = gmp_.quatProd(gmp_.quatExp(e0), Q0);
% T = Timed(end) / temp_s;
T = 8.33; % set the time duration of the executed motion

dt = Ts; % time step for numerical integration

% Execute DMP with new target orientation and time duration
[Time, Q_data, vRot_data, dvRot_data] = simulateGMPo(gmp_o, Q0, Qg, T, dt);
toc(t_start)


%% Groundtruth trajectory that should be produced
temp_s = Timed(end) / T; % temporal scaling 
Ks = gmp_o.getScaling(); % spatial scaling

n_data = size(Qd_data,2);

qLog_demo_data = zeros(3, n_data); % quatLog of the demo
qLog2_data = zeros(3, n_data); % quatLog of the groundtruth
% groundtruth orientation trajectory:
Q2_data = zeros(4, n_data); 
vRot2_data = zeros(3, n_data);
dvRot2_data = zeros(3, n_data);

for j=1:n_data
    qLog_demo_data(:,j) = GMPo.quat2q(Qd_data(:,j), Qd0);
    qLog2_data(:,j) = Ks*qLog_demo_data(:,j);
    Q2_data(:,j) = GMPo.q2quat(qLog2_data(:,j), Q0);
end

Time2 = Timed/temp_s;
dTime = diff(Time2);

for j=1:size(vRot2_data,2)-1
    vRot2_data(:,j) = gmp_.quatLogDiff(Q2_data(:,j+1),Q2_data(:,j)) / dTime(j);
end
vRot2_data(:,j) = zeros(3,1);

for i=1:3, dvRot2_data(i,:)=[diff(vRot2_data(i,:))./dTime 0]; end

% calc the quatLog of the DMP generated orientation
qLog_data = zeros(3, size(Q_data,2));
for j=1:size(qLog_data,2)
    qLog_data(:,j) = GMPo.quat2q(Q_data(:,j), Q0);
end


%% Plot results
line_width = 2.5;

figure('Position', [200 200 600 500]);
y_labels = {'$e_{q,x}$','$e_{q,y}$', '$e_{q,z}$'};
for i=1:3
   subplot(3,1,i);
   hold on;
   plot(Time, qLog_data(i,:), 'LineWidth', line_width, 'Color','blue');
   plot(Time2, qLog2_data(i,:), 'LineWidth', line_width, 'LineStyle',':', 'Color',[0.85 0.33 0.1 0.85]);
   plot(Timed, qLog_demo_data(i,:), 'LineWidth', line_width, 'LineStyle','--', 'Color',[0.93 0.69 0.13 0.7]);
   ylabel(y_labels{i}, 'interpreter','latex', 'fontsize',20);
   axis tight;
   if (i==1), legend({'$DMP$', '$ground-truth$', '$demo$'}, 'interpreter','latex', 'fontsize',16, 'Position',[0.7 0.78 0.27 0.15]); end
   if (i==1), title('Quaternion error: $e_q = log(Q * Q_0^{-1})$', 'interpreter','latex', 'fontsize',18); end
   if (i==3), xlabel('time [$s$]', 'interpreter','latex', 'fontsize',17); end
   hold off;
end

fig = figure;
fig.Position(3:4) = [659 606];
hold on;
plot3(qLog_data(1,:), qLog_data(2,:), qLog_data(3,:), 'LineWidth', line_width, 'LineStyle','-', 'Color','blue');
plot3(qLog2_data(1,:), qLog2_data(2,:), qLog2_data(3,:), 'LineWidth', line_width, 'LineStyle',':', 'Color',[0.85 0.33 0.1]);
plot3(qLog_demo_data(1,:), qLog_demo_data(2,:), qLog_demo_data(3,:), 'LineWidth', line_width, 'LineStyle','--', 'Color',[0.93 0.69 0.13]);
plot3(qLog2_data(1,1), qLog2_data(2,1), qLog2_data(3,1), 'LineWidth',4, 'Marker','o', 'MarkerSize',10, 'LineStyle','none', 'Color','green', 'DisplayName', '$p_0$');
plot3(qLog2_data(1,end), qLog2_data(2,end), qLog2_data(3,end), 'LineWidth',4, 'Marker','x', 'MarkerSize',10, 'LineStyle','none', 'Color','red', 'DisplayName', 'target');
plot3(qLog_demo_data(1,end), qLog_demo_data(2,end), qLog_demo_data(3,end), 'LineWidth',4, 'Marker','x', 'MarkerSize',10, 'LineStyle','none', 'Color','magenta', 'DisplayName', 'demo target');
legend({'$gmp$', '$ground-truth$', '$demo$'}, 'interpreter','latex', 'fontsize',17);
grid on;
view(-14.7, 58.8);
hold off;


figure;
Q_labels = {'$w$','$x$', '$y$', '$z$'};
Qd_labels = {'$w_d$','$x_d$', '$y_d$', '$z_d$'};
for i=1:4
   subplot(4,1,i);
   hold on;
   plot(Time, Q_data(i,:), 'LineWidth', line_width);
   plot(Time2, Q2_data(i,:), 'LineWidth', line_width, 'LineStyle',':');
   legend({Q_labels{i}, Qd_labels{i}}, 'interpreter','latex', 'fontsize',15);
   if (i==1), title('Unit Quaternion', 'interpreter','latex', 'fontsize',17); end
   if (i==4), xlabel('time [$s$]', 'interpreter','latex', 'fontsize',15); end
   hold off;
end

figure;
vRot_labels = {'$\omega_x$','$\omega_y$', '$\omega_z$'};
vRotd_labels = {'$\omega_{d,x}$','$\omega_{d,y}$', '$\omega_{d,z}$'};
for i=1:3
   subplot(3,1,i);
   hold on;
   plot(Time, vRot_data(i,:), 'LineWidth', line_width);
   plot(Time2, vRot2_data(i,:), 'LineWidth', line_width, 'LineStyle',':');
   legend({vRot_labels{i}, vRotd_labels{i}}, 'interpreter','latex', 'fontsize',15);
   if (i==1), title('Rotational Velocity', 'interpreter','latex', 'fontsize',17); end
   if (i==3), xlabel('time [$s$]', 'interpreter','latex', 'fontsize',15); end
   hold off;
end

figure;
dvRot_labels = {'$\dot{\omega}_x$','$\dot{\omega}_y$', '$\dot{\omega}_z$'};
dvRotd_labels = {'$\dot{\omega}_{d,x}$','$\dot{\omega}_{d,y}$', '$\dot{\omega}_{d,z}$'};
for i=1:3
   subplot(3,1,i);
   hold on;
   plot(Time, dvRot_data(i,:), 'LineWidth', line_width);
   plot(Time2, dvRot2_data(i,:), 'LineWidth', line_width, 'LineStyle',':');
   legend({dvRot_labels{i}, dvRotd_labels{i}}, 'interpreter','latex', 'fontsize',15);
   if (i==1), title('Rotational Acceleration', 'interpreter','latex', 'fontsize',17); end
   if (i==3), xlabel('time [$s$]', 'interpreter','latex', 'fontsize',15); end
   hold off;
end

%% ================================================
%% ================================================

function [Time, Q_data, rotVel_data, rotAccel_data] = simulateGMPo_in_Cart_space(gmp_o, Q0, Qg, T, dt)
%% Simulates a dmp encoding Cartesian orientation using unit quaternions.
% The DMP dynamics are formulated in the Cartesian orientation space.
% All quaternions are assumed to be unit, with format [w; x; y; z]
% @param[in] gmp: object of type @GMPo.
% @param[in] Q0: Initial orientation as a unit quaternion.
% @param[in] Qg: Target orientation as a unit quaternion.
% @param[in] T: Movement total time duration.
% @param[in] dt: Simulation timestep.
% @param[out] Time: 1 x N rowvector with simulation output timestamps.
% @param[out] Q_data: 4 x N matrix with simulation output positions.
% @param[out] rotVel_data: 3 x N matrix with simulation output rotational velocities.
% @param[out] rotAccel_data: 3 x N matrix with simulation output rotational accelerations.
%

%% data to log
Time = [];
Q_data = [];
rotVel_data = [];
rotAccel_data = [];

%% set initial values
t_end = T;
tau = t_end;

t = 0.0;

% phase variable, from 0 to 1
x = 0.0;
x_dot = 1/tau;
x_ddot = 0;

% Cartesian orientation trajectory
Q = Q0;
rotVel = zeros(3,1);
rotAccel = zeros(3,1);

gmp_o.setQ0(Q0);   % set initial orientation
gmp_o.setQg(Qg);   % set target orientation

% elaps_t = [];

%% simulate
while (true)

    %% data logging
    Time = [Time t];
    Q_data = [Q_data Q];
    rotVel_data = [rotVel_data rotVel];
    rotAccel_data = [rotAccel_data rotAccel];

%     t_start = tic;
    
    % Update the target if you want...
    % gmp.setQg(Qg);

    %% GMP simulation
    
    % get the scaled DMP learned orientation trajectory in the Cartesian
    % orientation space
    Qd = gmp_o.getQd(x);
    Vd = gmp_o.getVd(x, x_dot);
    Vd_dot = gmp_o.getVdDot(x, x_dot, x_ddot);
    % Or get them all together:
    %[Qd, Vd, Vd_dot] = gmp_o.getRefTraj(x, x_dot, x_ddot);

    external_signal = 0; % optionally add some external signal
    
    % Track it using a 2nd order dynamical system. 
    % DMP dynamics in the Cartesian orientation space:
    rotAccel = Vd_dot + 5*(Vd-rotVel) + 20*gmp_.quatLog(gmp_.quatDiff(Qd,Q));
    
    % if external_signal == 0 you can obviously set directly:
%     Q = Qd;
%     rotVel = Vd;
%     rotAccel = Vd_dot;

%     elaps_t = [elaps_t toc(t_start)*1000];

    %% Stopping criteria
    if (t>1.05*t_end)
        warning('Time limit reached... Stopping simulation!');
        break;
    end

    eo = gmp_.quatLog( gmp_.quatDiff(Qg,Q) );
    if (t>=t_end && norm(eo)<0.01)
        break;
    end

    %% Numerical integration
    t = t + dt;
    x = x + x_dot*dt;
    x_dot = x_dot + x_ddot*dt;
    Q = gmp_.quatProd( gmp_.quatExp(rotVel*dt), Q);
    rotVel = rotVel + rotAccel*dt;

end

% mean_elaps_t = mean(elaps_t)
% std_elaps_t = std(elaps_t,1)

end

function [Time, Q_data, rotVel_data, rotAccel_data] = simulateGMPo_in_log_space(gmp_o, Q0, Qg, T, dt)
%% Simulates a dmp encoding Cartesian orientation usning unit quaternions.
% The DMP dynamics are formulated in the quaterion-logarith space and the
% generated trajectory is converted to the Cartesian orientation space.
% All quaternions are assumed to be unit, with format [w; x; y; z]
% @param[in] gmp: object of type @GMPo.
% @param[in] Q0: Initial orientation as a unit quaternion.
% @param[in] Qg: Target orientation as a unit quaternion.
% @param[in] T: Movement total time duration.
% @param[in] dt: Simulation timestep.
% @param[out] Time: 1 x N rowvector with simulation output timestamps.
% @param[out] Q_data: 4 x N matrix with simulation output positions.
% @param[out] rotVel_data: 3 x N matrix with simulation output rotational velocities.
% @param[out] rotAccel_data: 3 x N matrix with simulation output rotational accelerations.
%

% data to log
Time = [];
Q_data = [];
rotVel_data = [];
rotAccel_data = [];

%% set initial values
t_end = T;
tau = t_end;

t = 0.0;

% phase variable, from 0 to 1
x = 0.0;
x_dot = 1/tau;
x_ddot = 0;

% Cartesian orientation trajectory
Q = Q0;
Q_prev = Q;
rotVel = zeros(3,1);
rotAccel = zeros(3,1);

% quaternion logarithm trajectory
y = gmp_o.quat2q(Q0, Q0);
dy = zeros(3,1);
ddy = zeros(3,1);

gmp_o.setQ0(Q0);   % set initial orientation
gmp_o.setQg(Qg);   % set target orientation

% elaps_t = [];

%% simulate
while (true)

    %% data logging
    Time = [Time t];
    Q_data = [Q_data Q];
    rotVel_data = [rotVel_data rotVel];
    rotAccel_data = [rotAccel_data rotAccel];
    
%     t_start = tic;

    %% DMP simulation
    
    % Deprecated, old DMP style evolution:
%     yc = zeros(3,1);
%     zc = zeros(3,1);
%     s = GMP_phase(x, x_dot, x_ddot);
%     gmp_o.update(s, y, z, yc, zc);
% 
%     dy = gmp_o.getYdot();
%     dz = gmp_o.getZdot();
%     rotAccel = gmp_o.getRotAccel(Q, yc_dot);

    % get the scaled DMP learned trajectory in the quatLog space.
    y_x = gmp_o.getYd(x);
    dy_x = gmp_o.getYdDot(x, x_dot);
    ddy_x = gmp_o.getYdDDot(x, x_dot, x_ddot);
    
    K = 20; % set the DMP stiffness in the quatLog space
    D = 5; % set the DMP damping in the quatLog space
    
    external_signal = 0; % optionally add some external signal
    
    % Track it using a 2nd order dynamical system.
    % DMP dynamics in the quatLog space:
    ddy = ddy_x + D*(dy_x - dy) + K*(y_x - y) + external_signal;
    
    % if external_signal == 0 you can obviously set directly:
%     y = y_x;
%     dy = dy_x;
%     ddy = ddy_x;

%     elaps_t = [elaps_t toc(t_start)*1000];
    
    %% Stopping criteria
    if (t>1.05*t_end)
        warning('Time limit reached... Stopping simulation!');
        break;
    end

    eo = gmp_.quatLog( gmp_.quatDiff(Qg,Q) );
    if (t>=t_end && norm(eo)<0.01)
        break;
    end

    %% Numerical integration
    t = t + dt;
    x = x + x_dot*dt;
    x_dot = x_dot + x_ddot*dt;
    y = y + dy*dt;
    dy = dy + ddy*dt;
    
    %% Convert trajectory from quatLog to the Cartesian orientation space
    Q = gmp_o.q2quat(y, Q0);
    % This is required to avoid discontinuities (since Q and -Q are the same orientation)
    if (dot(Q_prev,Q)<0)
        Q = -Q; 
    end
    Q_prev = Q;
    
    Q1 = gmp_o.getQ1(Q, Q0);
    rotVel = gmp_.qLogDot_to_rotVel(dy, Q1);
    rotAccel = gmp_.qLogDDot_to_rotAccel(ddy, rotVel, Q1);

end

% mean_elaps_t = mean(elaps_t)
% std_elaps_t = std(elaps_t,1)

end

