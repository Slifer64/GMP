%% Simulates a GMP
%  Loads a reference trajecory.
%  Trains a GMP based on the reference trajectory.
%  Plots and compares the results.

clc;
close all;
clear;

set_matlab_utils_path('../');

% fid = FileIO('data/train_data.bin', FileIO.in);
% Timed = fid.read('Timed');
% Pd_data = fid.read('Pd_data');
% dPd_data = fid.read('dPd_data');

Ts = 0.002;
[Timed, Pd_data, dPd_data, ddPd_data] = createData(0, 10, Ts, [0; 0; 0.8], [0.6; 0.6; 0.4], [0; 0; -0.8]);

% fid = FileIO('pos_data.bin', bitor(FileIO.out, FileIO.trunc) );
% fid.write('Timed',Timed);
% fid.write('Pd_data',Pd_data);
% fid.write('dPd_data',dPd_data);
% fid.write('ddPd_data',ddPd_data);
% fid.close();

% Data = struct('Time',Timed, 'Pos',Pd_data, 'Vel',dPd_data, 'Accel',ddPd_data);
% % save('pos_data.mat','Data');

% load('data/train_data.mat', 'Data');
% Timed = Data.Time;
% Pd_data = Data.Pos;
% dPd_data = Data.Vel;
% ddPd_data = zeros(size(dPd_data));
% dTime = diff(Timed);
% for i=1:size(ddPd_data,1)
%     ddPd_data(i,:) = [diff(dPd_data(i,:)) ./ dTime, 0];
% end

Ts = Timed(2) - Timed(1);
 
% %% Load training data
% load('data/train_data.mat', 'Data');
% 
% Timed = Data.Time;
% Pd_data = Data.Pos;
% dPd_data = Data.Vel;
% ddPd_data = Data.Accel;
% 
% Ts = Timed(2)-Timed(1);

read_from_file = false;

if (read_from_file)
    
    gmp = GMP_nDoF();
    gmp_.read(gmp, 'gmp_ndof.bin', 't1_')

else
    
    %% initialize and train GMP
    train_method = 'LS';
    N_kernels = 25;
    kernels_std_scaling = 1.5;
    n_dof = 3;
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

% gmp_.write(gmp, 'gmp_ndof.bin', 't1_');

%% DMP simulation
disp('GMP simulation...');
tic

spat_s = [2; 3; 1.1]; % spatial scale
temp_s = 1.2; % temporal scale
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

% Data.Time = Time;
% Data.Pos = P_data;
% Data.Vel = dP_data;
% Data.Accel = ddP_data;
% save('data/train_data.mat', 'Data');


function [Timed, yd_data, dyd_data, ddyd_data] = createData(t0, tf, Ts, p0, pf, pm)

    Timed = t0:Ts:tf;
    
    yd_data = get5thOrderPol(p0, pf, Timed);

    tm = [0.5 0.5 0.5]*Timed(end);
    sigma_m = [0.9 1.1 1.2];

    ym = zeros(size(yd_data));
    for i=1:size(ym,1)
       ym(i,:) = pm(i)*exp(-0.5*(Timed-tm(i)).^2/sigma_m(i)^2);
    end
%     ym = pm*exp(-0.5*(Timed-tm).^2/0.9^2);
    
    yd_data = yd_data + ym;
    dyd_data = zeros(size(yd_data));
    ddyd_data = zeros(size(yd_data));
    for i=1:3
        dyd_data(i,:) = [diff(yd_data(i,:)) 0]/Ts; 
        ddyd_data(i,:) = [diff(dyd_data(i,:)) 0]/Ts; 
    end

end
