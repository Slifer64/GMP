clc;
% close all;
clear;

import_io_lib();

% ------------------------------------------

fid = FileIO('data/orient_gmp_test_results.bin', FileIO.in);

Timed = fid.read('Timed');
Qd_data = fid.read('Qd_data');
% vRotd_data = fid.read('vRotd_data');
% dvRotd_data = fid.read('dvRotd_data');
Time = fid.read('Time');
Q_data = fid.read('Q_data');
vRot_data = fid.read('vRot_data');
dvRot_data = fid.read('dvRot_data');
Ks = fid.read('Ks');
temp_s = fid.read('temp_s');


%% Groundtruth trajectory that should be produced
n_data = size(Qd_data,2);

qLog_demo_data = zeros(3, n_data); % quatLog of the demo
qLog2_data = zeros(3, n_data); % quatLog of the groundtruth
% groundtruth orientation trajectory:
Q2_data = zeros(4, n_data); 
vRot2_data = zeros(3, n_data);
dvRot2_data = zeros(3, n_data);

Qd0 = Qd_data(:,1);
Q0 = Q_data(:,1);

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

