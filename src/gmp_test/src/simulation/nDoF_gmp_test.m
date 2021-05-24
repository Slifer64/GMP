%% Simulates a GMP

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

train_filename = 'data/pos_data.bin';

scale_type = 'rot_wb'; % {"prop", "rot_min", "rot_wb"}
wb_normal = [0; 0; 1]; % work-bench normal

spat_s = [2; 3; 1.1]; % spatial scale
temp_s = 1.2; % temporal scale

read_gmp_from_file = false;
write_gmp_to_file = false;

gmp_filename = 'data/gmp_pos.bin';
          

%% =============  Load train data  =============
fid = FileIO(train_filename, FileIO.in );
Timed = fid.read('Timed');
Pd_data = fid.read('Pd_data');
dPd_data = fid.read('dPd_data');
ddPd_data = fid.read('ddPd_data');
fid.close();

Ts = Timed(2) - Timed(1);


%% =============  Create/Train GMP  =============

n_dof = size(Pd_data, 1);

if (read_gmp_from_file)
    
    gmp = GMP_nDoF();
    gmp_.read(gmp, gmp_filename, '');

else
    
    %% initialize and train GMP
    gmp = GMP_nDoF(n_dof, N_kernels, kernels_std_scaling);
    tic
    offline_train_mse = gmp.train(train_method, Timed/Timed(end), Pd_data);
    offline_train_mse
    toc

    % traj_sc = TrajScale_Prop(n_dof);
    % traj_sc = TrajScale_Rot_min();
    traj_sc = TrajScale_Rot_wb();
    traj_sc.setWorkBenchNormal([0; 0; 1]);

    gmp.setScaleMethod(traj_sc);

end

% set scaling type
if ( strcmpi(scale_type, 'prop') ), traj_sc = TrajScale_Prop(n_dof);
elseif ( strcmpi(scale_type, 'rot_min') ), traj_sc = TrajScale_Rot_min();
elseif ( strcmpi(scale_type, 'rot_wb') )
    traj_sc = TrajScale_Rot_wb();
    traj_sc.setWorkBenchNormal(wb_normal);
else, error(['Unsupported scale type ''' scale_type '''...\n']);
end

gmp.setScaleMethod(traj_sc);
  

if (write_gmp_to_file), gmp_.write(gmp, gmp_filename, ''); end

%% DMP simulation
disp('GMP simulation...');
tic

P0d = Pd_data(:,1);
Pgd = Pd_data(:,end);
P0 = P0d;
Pg = spat_s.*(Pgd - P0) + P0;
T = Timed(end) / temp_s;
dt = Ts;

[Time, P_data, dP_data, ddP_data] = simulateGMP_nDoF(gmp, P0, Pg, T, dt);
toc

%% Reference trajectory (scaled)
Timed = Timed / temp_s;
Pd2_data = spat_s.*( Pd_data-P0 ) + P0;
dPd2_data = spat_s.*dPd_data*temp_s;
ddPd2_data = spat_s.*ddPd_data*temp_s^2;

%% Plot results
for i=1:3
    figure;
    subplot(3,1,1);
    hold on;
    plot(Time, P_data(i,:), 'LineWidth',2.0 , 'Color','blue');
    plot(Timed, Pd2_data(i,:), 'LineWidth',2.0, 'LineStyle',':', 'Color','magenta');
    ylabel('pos [$m$]', 'interpreter','latex', 'fontsize',15);
    title(['temporal scale: $' num2str(temp_s) '$     ,     spatial scale: $' num2str(spat_s(1)) ', ' num2str(spat_s(2)), ', ' num2str(spat_s(3)) '$'], 'interpreter','latex', 'fontsize',18);
    legend({'sim','$k_s$*demo'}, 'interpreter','latex', 'fontsize',15);

    axis tight;
    hold off;

    subplot(3,1,2);
    hold on;
    plot(Time, dP_data(i,:), 'LineWidth',2.0, 'Color','blue');
    plot(Timed, dPd2_data(i,:), 'LineWidth',2.0, 'LineStyle',':', 'Color','magenta');
    ylabel('vel [$m/s$]', 'interpreter','latex', 'fontsize',15);
    axis tight;
    hold off;

    subplot(3,1,3);
    hold on;
    plot(Time, ddP_data(i,:), 'LineWidth',2.0, 'Color','blue');
    plot(Timed, ddPd2_data(i,:), 'LineWidth',2.0, 'LineStyle',':', 'Color','magenta');
    ylabel('accel [$m/s^2$]', 'interpreter','latex', 'fontsize',15);
    axis tight;
    hold off;
end


figure;
hold on;
plot3(P_data(1,:), P_data(2,:), P_data(3,:), 'LineWidth',2, 'LineStyle','-', 'Color','magenta');
plot3(Pd2_data(1,:), Pd2_data(2,:), Pd2_data(3,:), 'LineWidth',2, 'LineStyle',':', 'Color','blue');
plot3(Pd_data(1,:), Pd_data(2,:), Pd_data(3,:), 'LineWidth',2, 'LineStyle',':', 'Color',[0 0.7 0]);
legend({'sim','$k_s$*demo','demo'}, 'interpreter','latex', 'fontsize',15);
xlabel('$x$', 'interpreter','latex', 'fontsize',15);
ylabel('$y$', 'interpreter','latex', 'fontsize',15);
zlabel('$z$', 'interpreter','latex', 'fontsize',15);
hold off;


%% ============================================================
%% ============================================================


function [Time, Y_data, dY_data, ddY_data] = simulateGMP_nDoF(gmp, y0, g, T, dt)
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
Dim = gmp.numOfDoFs();
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
s = GMP_phase(x, x_dot, 0);

iters = 0;
Time = [];
Y_data = [];
dY_data = [];
ddY_data = [];
x_data = [];

gmp.setY0(y0);
gmp.setGoal(g);

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

    %% Stopping criteria
    if (t>=1.1*t_end) % && norm(y-g)<5e-3 && norm(dy)<5e-3)
        break;
    end

    %% Numerical integration
    iters = iters + 1;
    t = t + dt;
    x = x + x_dot*dt;
    y = y + dy*dt;
    z = z + dz*dt;

    s.x = x;
    
end


end
