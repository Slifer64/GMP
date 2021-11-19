clc;
close all;
clear;

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
Ks = fid.read('Ks');
temp_s = fid.read('temp_s');
% % Or just:
% data = fid.readAll();
% Timed = data.Timed 
% ...

P0d = Pd_data(:,1);
P0 = P_data(:,1);

%% This is the groundtruth trajectory that should be produced
Timed2 = Timed / temp_s;
Pd2_data = Ks*( Pd_data-P0d ) + P0;
dPd2_data = Ks*dPd_data*temp_s;
ddPd2_data = Ks*ddPd_data*temp_s^2;

Pgd = Pd_data(:,end);
Pg = Pd2_data(:,end);

%% Plot results
for i=1:3
    figure;
    subplot(3,1,1);
    hold on;
    plot(Time, P_data(i,:), 'LineWidth',2.0 , 'Color','blue', 'DisplayName', 'DMP');
    plot(Timed2, Pd2_data(i,:), 'LineWidth',2.0, 'LineStyle',':', 'Color','magenta', 'DisplayName','ground-truth');
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
