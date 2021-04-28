%% Simulates a GMP as a DMP
%  Loads a reference trajecory.
%  Trains a GMP based on the reference trajectory.
%  Plots and compares the results.
function test_nDoF_gmp()

set_matlab_utils_path();


fid = FileIO('data/train_data.bin', FileIO.in);
% fid.printHeader();

Timed = fid.read('Timed');
Pd_data = fid.read('Pd_data');
dPd_data = fid.read('dPd_data');
% Qd_data = fid.read('Qd_data');
% vRotd_data = fid.read('vRotd_data');
% joint_pos_data = fid.read('joint_pos_data');

ddPd_data = zeros(size(dPd_data));
dTime = diff(Timed);
for i=1:size(ddPd_data,1)
    ddPd_data(i,:) = [diff(dPd_data(i,:)) ./ dTime, 0];
end

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

%% initialize and train GMP
train_method = 'LS';
N_kernels = 30;
kernels_std_scaling = 1;
n_dof = 3;
gmp = GMP_nDoF(n_dof, N_kernels, 30, 100, kernels_std_scaling);
tic
offline_train_mse = gmp.train(train_method, Timed, Pd_data);
offline_train_mse
toc

% gmp.exportToFile('data/gmp_nDoF_model.bin');
% gmp = GMP_nDoF.importFromFile('data/gmp_nDoF_model.bin');

%% DMP simulation
disp('GMP simulation...');
tic

spat_s = 1; % spatial scale
temp_s = 1; % temporal scale
P0 = Pd_data(:,1);
Pgd = Pd_data(:,end);
Pg = P0 + spat_s*(Pgd - P0);
T = Timed(end) / temp_s;
dt = Ts;
[Time, P_data, dP_data, ddP_data] = simulateGMP_nDoF(gmp, P0, Pg, T, dt);
toc


%% Reference trajectory (scaled)
Timed = Timed / temp_s;
Pd_data = spat_s*( Pd_data-P0 ) + P0;
dPd_data = spat_s*dPd_data*temp_s;
ddPd_data = spat_s*ddPd_data*temp_s^2;

%% Plot results
for i=1:3
    figure;
    subplot(3,1,1);
    hold on;
    plot(Time, P_data(i,:), 'LineWidth',2.0 , 'Color','blue');
    plot(Timed, Pd_data(i,:), 'LineWidth',2.0, 'LineStyle',':', 'Color','magenta');
    ylabel('pos [$m$]', 'interpreter','latex', 'fontsize',15);
    title(['temporal scale: $' num2str(temp_s) '$     ,     spatial scale: $' num2str(spat_s) '$'], 'interpreter','latex', 'fontsize',18);
    legend({'sim','demo'}, 'interpreter','latex', 'fontsize',15);

    axis tight;
    hold off;

    subplot(3,1,2);
    hold on;
    plot(Time, dP_data(i,:), 'LineWidth',2.0, 'Color','blue');
    plot(Timed, dPd_data(i,:), 'LineWidth',2.0, 'LineStyle',':', 'Color','magenta');
    ylabel('vel [$m/s$]', 'interpreter','latex', 'fontsize',15);
    axis tight;
    hold off;

    subplot(3,1,3);
    hold on;
    plot(Time, ddP_data(i,:), 'LineWidth',2.0, 'Color','blue');
    plot(Timed, ddPd_data(i,:), 'LineWidth',2.0, 'LineStyle',':', 'Color','magenta');
    ylabel('accel [$m/s^2$]', 'interpreter','latex', 'fontsize',15);
    axis tight;
    hold off;
end


figure;
hold on;
plot3(P_data(1,:), P_data(2,:), P_data(3,:), 'LineWidth',2, 'LineStyle','-', 'Color','magenta');
plot3(Pd_data(1,:), Pd_data(2,:), Pd_data(3,:), 'LineWidth',2, 'LineStyle',':', 'Color','blue');
hold off;

end


