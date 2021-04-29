clc;
close all;
clear;

addpath('../../matlab/lib/io_lib/');
import_io_lib();

% ------------------------------------------

fid = FileIO('data/nDoF_gmp_test_results.bin', FileIO.in);
Timed = fid.read('Timed');
Pd_data = fid.read('Pd_data');
dPd_data = fid.read('dPd_data');
ddPd_data = fid.read('ddPd_data');
Time = fid.read('Time');
P_data = fid.read('P_data');
dP_data = fid.read('dP_data');
ddP_data = fid.read('ddP_data');
spat_s = fid.read('spat_s');
temp_s = fid.read('temp_s');

%% Reference trajectory (scaled)
P0 = Pd_data(:,1);
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
