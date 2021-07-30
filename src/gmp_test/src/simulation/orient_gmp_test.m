% function orient_gmp_test()

clc;
close all;
clear;

%% =============  includes...  =============
addpath('../../../matlab/lib/gmp_lib/');
import_gmp_lib();

addpath('../../../matlab/lib/io_lib/');
import_io_lib();

%% =============  Load params  =============
train_method = 'LS';
N_kernels = 25;
kernels_std_scaling = 1.5;

train_filename = 'data/orient_data.bin';

scale_type = 'rot_wb'; % {"prop", "rot_min", "rot_wb"}
wb_normal = [0; 0; 1]; % work-bench normal

spat_s = [2; 3; 1.1]; % spatial scale
temp_s = 1.2; % temporal scale

read_gmp_from_file = false;
write_gmp_to_file = false;

gmp_filename = 'data/gmp_orient.bin';


%% Load training data
fid = FileIO(train_filename, FileIO.in );
Timed = fid.read('Timed');
Qd_data = fid.read('Qd_data');
vRotd_data = fid.read('vRotd_data');
dvRotd_data = fid.read('dvRotd_data');
fid.close();

Ts = Timed(2)-Timed(1);

simulateGMPo = @simulateGMPo_in_Cart_space; % simulateGMPo_in_'log/quat/Cart'_space

read_from_file = true;

if (read_from_file)
    
    gmp_o = GMPo();
    gmp_.read(gmp_o, gmp_filename, '');

else
    
    %% initialize and train GMP
    train_method = 'LS';
    N_kernels = 30;
    kernels_std_scaling = 1;
    gmp_o = GMPo(N_kernels, kernels_std_scaling);
    tic
    offline_train_mse = gmp_o.train(train_method, Timed/Timed(end), Qd_data);
    offline_train_mse
    toc

    % traj_sc = TrajScale_Prop(n_dof);
    % traj_sc = TrajScale_Rot_min();
    traj_sc = TrajScale_Rot_wb();
    traj_sc.setWorkBenchNormal([0; 0; 1]);
    gmp_o.setScaleMethod(traj_sc);
    
end

% gmp_.write(gmp_o, 'gmp_o.bin', 'o_');

%% DMP simulation
disp('GMP simulation...');
tic
Qd0 = Qd_data(:,1);
Q0 = Qd0; %gmp_.quatProd(rotm2quat(rotz(57))',Qd0);
Qgd = Qd_data(:,end);
ks = diag([1.2 1.3 1.1]);
kt = 2;
e0 = ks*gmp_.quatLog(gmp_.quatDiff(Qgd,Qd0));
Qg = gmp_.quatProd(gmp_.quatExp(e0), Q0);
T = Timed(end) / kt;
dt = Ts;

% gmp_o.setScaleMethod(TrajScale.PROP_SCALE);

[Time, Q_data, vRot_data, dvRot_data] = simulateGMPo(gmp_o, Q0, Qg, T, dt);
toc


%% Plot results

Timed = Timed/kt;

T_sc = gmp_o.getScaling();

Pqd_data0 = zeros(3, size(Qd_data,2));
Pqd_data = zeros(3, size(Qd_data,2));
for j=1:size(Pqd_data,2)
    Pqd_data0(:,j) = GMPo.quat2q(Qd_data(:,j), Qd0);
    Pqd_data(:,j) = T_sc*Pqd_data0(:,j);
    Qd_data(:,j) = GMPo.q2quat(Pqd_data(:,j), Q0);
end

for j=1:size(vRotd_data,2)-1
    vRotd_data(:,j) = gmp_.quatLog(gmp_.quatDiff(Qd_data(:,j+1),Qd_data(:,j)))/Ts;
end
vRotd_data(:,j) = zeros(3,1);
for i=1:3, dvRotd_data(i,:)=[diff(vRotd_data(i,:)) 0]/Ts; end


Pq_data = zeros(3, size(Q_data,2));
for j=1:size(Pq_data,2)
    Pq_data(:,j) = GMPo.quat2q(Q_data(:,j), Q0);
end

line_width = 2.5;

figure('Position', [200 200 600 500]);
y_labels = {'$e_{q,x}$','$e_{q,y}$', '$e_{q,z}$'};
for i=1:3
   subplot(3,1,i);
   hold on;
   plot(Time, Pq_data(i,:), 'LineWidth', line_width, 'Color','blue');
   plot(Timed, Pqd_data(i,:), 'LineWidth', line_width, 'LineStyle',':', 'Color',[0.85 0.33 0.1 0.85]);
   plot(Timed, Pqd_data0(i,:), 'LineWidth', line_width, 'LineStyle','--', 'Color',[0.93 0.69 0.13 0.7]);
   ylabel(y_labels{i}, 'interpreter','latex', 'fontsize',20);
   axis tight;
   if (i==1), legend({'$gmp$', '$k_s * demo$', '$demo$'}, 'interpreter','latex', 'fontsize',16, 'Position',[0.7 0.78 0.27 0.15]); end
   if (i==1), title('Quaternion error: $e_q = log(Q * Q_0^{-1})$', 'interpreter','latex', 'fontsize',18); end
   if (i==3), xlabel('time [$s$]', 'interpreter','latex', 'fontsize',17); end
   hold off;
end

figure;
hold on;
plot3(Pq_data(1,:), Pq_data(2,:), Pq_data(3,:), 'LineWidth', line_width, 'LineStyle','-', 'Color','blue');
plot3(Pqd_data(1,:), Pqd_data(2,:), Pqd_data(3,:), 'LineWidth', line_width, 'LineStyle',':', 'Color',[0.85 0.33 0.1]);
plot3(Pqd_data0(1,:), Pqd_data0(2,:), Pqd_data0(3,:), 'LineWidth', line_width, 'LineStyle','--', 'Color',[0.93 0.69 0.13]);
legend({'$gmp$', '$k_s * demo$', '$demo$'}, 'interpreter','latex', 'fontsize',17);
hold off;


figure;
Q_labels = {'$\eta$','$\epsilon_1$', '$\epsilon_2$', '$\epsilon_3$'};
Qd_labels = {'$\eta_d$','$\epsilon_{d,1}$', '$\epsilon_{d,2}$', '$\epsilon_{d,3}$'};
for i=1:4
   subplot(4,1,i);
   hold on;
   plot(Time, Q_data(i,:), 'LineWidth', line_width);
   plot(Timed, Qd_data(i,:), 'LineWidth', line_width, 'LineStyle',':');
   legend({Q_labels{i}, Qd_labels{i}}, 'interpreter','latex', 'fontsize',15);
   if (i==1), title('Unit Quaternion', 'interpreter','latex', 'fontsize',17); end
   if (i==4), xlabel('time [$s$]', 'interpreter','latex', 'fontsize',15); end
   hold off;
end

return

figure;
vRot_labels = {'$\omega_x$','$\omega_y$', '$\omega_z$'};
vRotd_labels = {'$\omega_{d,x}$','$\omega_{d,y}$', '$\omega_{d,z}$'};
for i=1:3
   subplot(3,1,i);
   hold on;
   plot(Time, vRot_data(i,:), 'LineWidth', line_width);
   plot(Timed, vRotd_data(i,:), 'LineWidth', line_width, 'LineStyle',':');
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
   plot(Timed, dvRotd_data(i,:), 'LineWidth', line_width, 'LineStyle',':');
   legend({dvRot_labels{i}, dvRotd_labels{i}}, 'interpreter','latex', 'fontsize',15);
   if (i==1), title('Rotational Acceleration', 'interpreter','latex', 'fontsize',17); end
   if (i==3), xlabel('time [$s$]', 'interpreter','latex', 'fontsize',15); end
   hold off;
end

%% ================================================
%% ================================================

function [Time, Q_data, rotVel_data, rotAccel_data] = simulateGMPo_in_Cart_space(gmp_o, Q0, Qg, T, dt)
%% Simulates a dmp encoding Cartesian orientation usning unit quaternions.


%% set initial values
t_end = T;
tau = t_end;

Time = [];
Q_data = [];
rotVel_data = [];
rotAccel_data = [];

t = 0.0;
x = 0.0;
x_dot = 1/tau;
x_ddot = 0;
Q = Q0;
rotVel = zeros(3,1);
rotAccel = zeros(3,1);

gmp_o.setQ0(Q0);   % set initial orientation
gmp_o.setQg(Qg);   % set target orientation

elaps_t = [];

%% simulate
while (true)

    %% data logging
    Time = [Time t];
    Q_data = [Q_data Q];
    rotVel_data = [rotVel_data rotVel];  
    rotAccel_data = [rotAccel_data rotAccel];

    tic;
    
    %% GMP simulation
%     Qd = gmp_o.getQd(x);
%     Vd = gmp_o.getVd(x, x_dot);
%     Vd_dot = gmp_o.getVdDot(x, x_dot, x_ddot);
    [Qd, Vd, Vd_dot] = gmp_o.getRefTraj(x, x_dot, x_ddot);
    
    rotAccel = Vd_dot + 5*(Vd-rotVel) + 20*gmp_.quatLog(gmp_.quatDiff(Qd,Q));
    
    elaps_t = [elaps_t toc()*1000];
    
    %% Stopping criteria   
    if (t>1.5*t_end)
        warning('Time limit reached... Stopping simulation!');
        break;
    end
    
    eo = gmp_.quatLog(gmp_.quatProd(Qg, gmp_.quatInv(Q)));
    if (t>=t_end && norm(eo)<0.02)
        break;
    end

    %% Numerical integration
    t = t + dt;
    x = x + x_dot*dt;
    Q = gmp_.quatProd( gmp_.quatExp(rotVel*dt), Q);
    rotVel = rotVel + rotAccel*dt;  
    
end

mean_elaps_t = mean(elaps_t)
std_elaps_t = std(elaps_t,1)

end

function [Time, Q_data, rotVel_data, rotAccel_data] = simulateGMPo_in_log_space(gmp_o, Q0, Qg, T, dt)
%% Simulates a dmp encoding Cartesian orientation usning unit quaternions.


%% set initial values
t_end = T;
tau = t_end;

Time = [];
Q_data = [];
rotVel_data = [];
rotAccel_data = [];

t = 0.0;
x = 0.0;
x_dot = 1/tau;
x_ddot = 0;
Q = Q0;
Q_prev = Q;
rotVel = zeros(3,1);
rotAccel = zeros(3,1);
q = gmp_o.quat2q(Q0, Q0);
qdot = zeros(3,1);
dy = zeros(3,1);
dz = zeros(3,1);

gmp_o.setQ0(Q0);
gmp_o.setQg(Qg);
y = gmp_o.getY(Q);
z = gmp_o.getZ(rotVel, Q);
g = gmp_o.quat2q(Qg, Q0);

elaps_t = [];

%% simulate
while (true)

    %% data logging
    Time = [Time t];
    Q_data = [Q_data Q];
    rotVel_data = [rotVel_data rotVel];  
    rotAccel_data = [rotAccel_data rotAccel];
    
    yc_dot = 0;
    
    tic

    %% DMP simulation
    yc = zeros(3,1);
    zc = zeros(3,1);
    s = GMP_phase(x, x_dot, x_ddot);
    gmp_o.update(s, y, z, yc, zc);

    dy = gmp_o.getYdot();
    dz = gmp_o.getZdot();
    rotAccel = gmp_o.getRotAccel(Q, yc_dot);
    
    elaps_t = [elaps_t toc()*1000];

    %% Stopping criteria   
    if (t>1.5*t_end)
        warning('Time limit reached... Stopping simulation!');
        break;
    end
    
    eo = quatLog(quatProd(Qg, quatInv(Q)));
    if (t>=t_end && norm(eo)<0.02)
        break;
    end

    %% Numerical integration
    t = t + dt;
    x = x + x_dot*dt;
    y = y + dy*dt;
    z = z + dz*dt;
    
    q = y;
    dy = z;
    qdot = dy;
    
    Q_prev = Q;
    Q = gmp_o.q2quat(q, Q0);
    if (Q_prev'*Q<0), Q = -Q; end
    
    Q1 = gmp_o.getQ1(Q, Q0);
    rotVel = gmp_o.qLogDot_to_rotVel(qdot, Q1);
    
end

mean_elaps_t = mean(elaps_t)
std_elaps_t = std(elaps_t,1)

end

function [Time, Q_data, rotVel_data, rotAccel_data] = simulateGMPo_in_quat_space(gmp_o, Q0, Qg, T, dt)
%% Simulates a dmp encoding Cartesian orientation usning unit quaternions.


%% set initial values
t_end = T;
tau = t_end;

Time = [];
Q_data = [];
rotVel_data = [];
rotAccel_data = [];

t = 0.0;
x = 0.0;
x_dot = 1/tau;
x_ddot = 0;
Q = Q0;
rotVel = zeros(3,1);
rotAccel = zeros(3,1);

gmp_o.setQ0(Q0);   % set initial orientation
gmp_o.setQg(Qg);   % set target orientation

elaps_t = [];

%% simulate
while (true)

    %% data logging
    Time = [Time t];
    Q_data = [Q_data Q];
    rotVel_data = [rotVel_data rotVel];  
    rotAccel_data = [rotAccel_data rotAccel];

    tic;
    
    %% GMP simulation
    yc = 0; % optional coupling for 'y' state
    zc = 0; % optional coupling for 'z' state
    yc_dot = 0; % derivative of coupling for 'y' state
    s = GMP_phase(x,x_dot,x_ddot);
    rotAccel = gmp_o.calcRotAccel(s, Q, rotVel, yc, zc, yc_dot);

    elaps_t = [elaps_t toc()*1000];
    
    %% Stopping criteria   
    if (t>1.5*t_end)
        warning('Time limit reached... Stopping simulation!');
        break;
    end
    
    eo = quatLog(quatProd(Qg, quatInv(Q)));
    if (t>=t_end && norm(eo)<0.02)
        break;
    end

    
    %% Numerical integration
    t = t + dt;
    x = x + x_dot*dt;
    Q = quatProd( quatExp(rotVel*dt), Q);
    rotVel = rotVel + rotAccel*dt;
    
end

mean_elaps_t = mean(elaps_t)
std_elaps_t = std(elaps_t,1)

end




