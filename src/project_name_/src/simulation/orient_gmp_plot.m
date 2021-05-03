clc;
close all;
clear;

addpath('../../matlab/lib/io_lib/');
import_io_lib();

% ------------------------------------------

fid = FileIO('data/orient_gmp_test_results.bin', FileIO.in);

Timed = fid.read('Timed');
Qd_data = fid.read('Qd_data');
vRotd_data = fid.read('vRotd_data');
dvRotd_data = fid.read('dvRotd_data');
Time = fid.read('Time');
Q_data = fid.read('Q_data');
vRot_data = fid.read('vRot_data');
dvRot_data = fid.read('dvRot_data');
T_sc = fid.read('T_sc');
kt = fid.read('temp_s');
    

%% Plot results

Timed = Timed/kt;

Qd0 = Qd_data(:,1);
Q0 = Q_data(:,1);

Pqd_data0 = zeros(3, size(Qd_data,2));
Pqd_data = zeros(3, size(Qd_data,2));
for j=1:size(Pqd_data,2)
    Pqd_data0(:,j) = GMPo.quat2q(Qd_data(:,j), Qd0);
    Pqd_data(:,j) = T_sc*Pqd_data0(:,j);
    Qd_data(:,j) = GMPo.q2quat(Pqd_data(:,j), Q0);
end

% for j=1:size(vRotd_data,2)-1
%     vRotd_data(:,j) = math_.quatLog(math_.quatDiff(Qd_data(:,j+1),Qd_data(:,j)))/Ts;
% end
% vRotd_data(:,j) = zeros(3,1);
% for i=1:3, dvRotd_data(i,:)=[diff(vRotd_data(i,:)) 0]/Ts; end


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

