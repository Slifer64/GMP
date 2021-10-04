clc;
close all;
clear;

%% =============  includes...  =============
import_gmp_lib();
import_io_lib();

%% =============  Load params  =============
train_method = 'LS';
N_kernels = 25;
kernels_std_scaling = 1.5;

train_filename = 'data/orient_data.bin';

scale_type = 'prop'; % {"prop", "rot_min", "rot_wb"}
wb_normal = [0; 0; 1]; % work-bench normal

spat_s = [1 1.2 1]; % spatial scale
temp_s = 1.5; % temporal scale


%% Load training data
fid = FileIO(train_filename, FileIO.in );
Timed = fid.read('Timed');
Qd_data = fid.read('Qd_data');
vRotd_data = fid.read('vRotd_data');
dvRotd_data = fid.read('dvRotd_data');
fid.close();

Ts = Timed(2)-Timed(1);

simulateGMPo = @simulateGMPo_in_Cart_space; % simulateGMPo_in_'log/quat/Cart'_space


    
%% initialize and train GMP
gmp_o = GMPo(N_kernels, kernels_std_scaling);
tic
offline_train_mse = gmp_o.train(train_method, Timed/Timed(end), Qd_data);
offline_train_mse
toc

if ( strcmpi(scale_type, 'prop') )
   traj_sc = TrajScale_Prop(3);
elseif ( strcmpi(scale_type, 'rot_min') )
    traj_sc = TrajScale_Rot_min();
elseif ( strcmpi(scale_type, 'rot_wb') )
    traj_sc = TrajScale_Rot_wb();
    traj_sc.setWorkBenchNormal(wb_normal);
else
    error(['Unsupported scale method ''' scale_type '''...']);
end

gmp_o.setScaleMethod(traj_sc);


%% DMP simulation
disp('GMP simulation...');
tic
Qd0 = Qd_data(:,1);
Q0 = Qd0; %gmp_.quatProd(rotm2quat(rotz(57))',Qd0);
Qgd = Qd_data(:,end);
ks = diag(spat_s);
kt = temp_s;
e0 = ks*gmp_.quatLogDiff(Qgd,Qd0);
Qg = gmp_.quatProd(gmp_.quatExp(e0), Q0);
% Qg = Q0 + [0.05; 0.08; 0.08; 0.07];
% Qg = Qg/norm(Qg);
T = Timed(end) / kt;
dt = Ts;

% gmp_o.setScaleMethod(TrajScale.PROP_SCALE);

%%  Optimization
x_data = 0:0.01:1;
theta_bar = 2.8;
qdot_bar = 2;
qddot_bar = 4;

% --------- optimize position ------------
tic
gmp_o2 = gmp_o.deepCopy();
[eflag, exit_msg] = optimizeGMPo(gmp_o2, Q0, Qg, T, x_data, theta_bar, qdot_bar, qddot_bar, 0, 1);
elaps_t = toc;
fprintf('===> Velocity Optimization finished:\n');
fprintf('Elaps time: %f ms\n',elaps_t*1000);
fprintf('Exit msg: %s\n',exit_msg);

[Time, Q_data, vRot_data, dvRot_data] = simulateGMPo(gmp_o2, Q0, Qg, T, dt);
[q_data, dq_data, ddq_data] = cartesian_to_quatLog(Q_data, vRot_data, dvRot_data);
q_data_norm = mat_col_norm(q_data);
dq_data_norm = mat_col_norm(dq_data);
ddq_data_norm = mat_col_norm(ddq_data);

% --------- optimize velocity ------------
tic
gmp_o2 = gmp_o.deepCopy();
[eflag, exit_msg] = optimizeGMPo(gmp_o2, Q0, Qg, T, x_data, theta_bar, qdot_bar, qddot_bar, 1, 0);
elaps_t = toc;
fprintf('===> Position Optimization finished:\n');
fprintf('Elaps time: %f ms\n',elaps_t*1000);
fprintf('Exit msg: %s\n',exit_msg);

[Time2, Q_data2, vRot_data2, dvRot_data2] = simulateGMPo(gmp_o2, Q0, Qg, T, dt);

% --------- prop ------------
[Time3, Q_data3, vRot_data3, dvRot_data3] = simulateGMPo(gmp_o, Q0, Qg, T, dt);
[qd_data, dqd_data, ddqd_data] = cartesian_to_quatLog(Q_data3, vRot_data3, dvRot_data3);
qd_data_norm = mat_col_norm(qd_data);
dqd_data_norm = mat_col_norm(dqd_data);
ddqd_data_norm = mat_col_norm(ddqd_data);
[q_data2, dq_data2, ddq_data2] = cartesian_to_quatLog(Q_data2, vRot_data2, dvRot_data2);
q_data_norm2 = mat_col_norm(q_data2);
dq_data_norm2 = mat_col_norm(dq_data2);
ddq_data_norm2 = mat_col_norm(ddq_data2);

% fid = FileIO(train_filename, bitor(FileIO.out,FileIO.trunc) );
% fid.write('Timed', Time);
% fid.write('Qd_data', Q_data);
% fid.write('vRotd_data', vRot_data);
% fid.write('dvRotd_data', dvRot_data);
% fid.close();
% return

%% Plot results
line_width = 2.5;
ax_fontsize = 15;
label_font = 18;

fig = figure; 
fig.Position(3:4) = [670 866];
ax = subplot(5,1,1); hold on;
plot(Time, q_data_norm, 'LineWidth', line_width, 'Color','green');
plot(Time, q_data_norm2, 'LineWidth', line_width, 'Color',[0.85 0.33 0.1]);
plot(Time, qd_data_norm, 'LineWidth', line_width, 'LineStyle',':', 'Color','blue');
ylabel('$||\eta||$ $[rad]$', 'interpreter','latex', 'fontsize',label_font);
axis tight;
plot(ax.XLim, [theta_bar theta_bar], 'LineWidth', line_width, 'LineStyle','--', 'Color','red');
ax.FontSize = ax_fontsize;
% -------------------------------------------
ax = subplot(5,1,2); hold on;
plot(Time, dq_data_norm, 'LineWidth', line_width, 'Color','green');
plot(Time, dq_data_norm2, 'LineWidth', line_width, 'Color',[0.85 0.33 0.1]);
plot(Time, dqd_data_norm, 'LineWidth', line_width, 'LineStyle',':', 'Color','blue');
ylabel('$||\dot{\eta}||$ $[rad/s]$', 'interpreter','latex', 'fontsize',label_font);
axis tight;
plot(ax.XLim, [qdot_bar qdot_bar], 'LineWidth', line_width, 'LineStyle','--', 'Color','red');
ax.FontSize = ax_fontsize;
% -------------------------------------------
ax = subplot(5,1,3); hold on;
plot(Time, ddq_data_norm, 'LineWidth', line_width, 'Color','green');
plot(Time, ddq_data_norm2, 'LineWidth', line_width, 'Color',[0.85 0.33 0.1]);
plot(Time, ddqd_data_norm, 'LineWidth', line_width, 'LineStyle',':', 'Color','blue');
ylabel('$||\ddot{\eta}||$ $[rad/s^2]$', 'interpreter','latex', 'fontsize',label_font);
% xlabel('time $[s]$', 'interpreter','latex', 'fontsize',label_font);
axis tight;
plot(ax.XLim, [qddot_bar qddot_bar], 'LineWidth', line_width, 'LineStyle','--', 'Color','red');
ax.FontSize = ax_fontsize;

yg = gmp_.quatLogDiff(Qg,Q0);
ygd = gmp_.quatLogDiff(Qgd,Q0);

% figure;
ax = subplot(5,1,[4 5]); hold on;
plot3(0, 0, 0, 'LineWidth', 4, 'LineStyle','none', 'Color',[0 0.7 0],'Marker','o', 'MarkerSize',12, 'HandleVisibility','off');
plot3(yg(1), yg(2), yg(3), 'LineWidth', 4, 'LineStyle','none', 'Color','red','Marker','x', 'MarkerSize',12, 'HandleVisibility','off');
plot3(ygd(1), ygd(2), ygd(3), 'LineWidth', 4, 'LineStyle','none', 'Color','magenta','Marker','x', 'MarkerSize',12, 'HandleVisibility','off');
plot3(q_data(1,:), q_data(2,:), q_data(3,:), 'LineWidth', line_width, 'LineStyle','-', 'Color','green');
plot3(q_data2(1,:), q_data2(2,:), q_data2(3,:), 'LineWidth', line_width, 'LineStyle','-', 'Color',[0.85 0.33 0.1]);
plot3(qd_data(1,:), qd_data(2,:), qd_data(3,:), 'LineWidth', line_width, 'LineStyle',':', 'Color','blue');
% plot3(qd_data0(1,:), qd_data0(2,:), qd_data0(3,:), 'LineWidth', line_width, 'LineStyle','--', 'Color',[0.93 0.69 0.13]);
legend({'$opt-vel$', '$opt-vel$', '$prop$'}, 'interpreter','latex', 'fontsize',17, 'Orientation','horizontal', 'Position',[0.1925 0.9410 0.6627 0.0366]);
grid on;
view(-70, 14.8);
ax.FontSize = ax_fontsize;
hold off;



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
    if (t>1.05*t_end)
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
    
    Q = Qd;
    rotVel = Vd;
    rotAccel = Vd_dot;
    
end

% mean_elaps_t = mean(elaps_t)
% std_elaps_t = std(elaps_t,1)

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
    if (t>1.05*t_end)
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
    if (t>1.05*t_end)
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

%% ===============================================
function [q_data, qdot_data, qddot_data] = cartesian_to_quatLog(Q_data, vRot_data, dvRot_data)
    
    Q0 = Q_data(:,1);
    n = size(Q_data,2);
    
    q_data = zeros(3,n);
    qdot_data = zeros(3,n);
    qddot_data = zeros(3,n);
    
    for j=1:n
        Q1 = gmp_.quatDiff(Q_data(:,j),Q0);
        q_data(:,j) = gmp_.quatLog(Q1);
        qdot_data(:,j) = gmp_.rotVel_to_qLogDot(vRot_data(:,j), Q1);
        qddot_data(:,j) = gmp_.rotAccel_to_qLogDDot(dvRot_data(:,j), vRot_data(:,j), Q1);
    end

end

function data_norm = mat_col_norm(data)

    data_norm = zeros(1, size(data,2));
    for i=1:length(data_norm)
        data_norm(i) = norm(data(:,i));
    end

end

