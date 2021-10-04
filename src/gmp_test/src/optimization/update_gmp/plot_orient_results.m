clc;
close all;
clear;

import_io_lib();

%% Load data
fid = FileIO('data/orient_results.bin', FileIO.in);
Timed           = fid.read('Timed');
q_data          = fid.read('q_data');
vRot_data       = fid.read('vRot_data');
dvRot_data      = fid.read('dvRot_data');

q_new_data      = fid.read('q_new_data');
vRot_new_data   = fid.read('vRot_new_data');
dvRot_new_data  = fid.read('dvRot_new_data');

qd_new_data     = fid.read('qd_new_data');
vRotd_new_data  = fid.read('vRotd_new_data');
dvRotd_new_data = fid.read('dvRotd_new_data');

up_ind          = fid.read('up_ind');


figure;
for i=1:3
   subplot(3,1,i); hold on;
   plot(Timed, qd_new_data(i,:), 'LineWidth',2, 'Color','blue', 'LineStyle','-');
   plot(Timed, q_data(i,:), 'LineWidth',2, 'Color',[0.85 0.33 0.1], 'LineStyle','--');
   plot(Timed, q_new_data(i,:), 'LineWidth',2, 'Color','magenta', 'LineStyle',':');
   plot(Timed(up_ind), qd_new_data(i,up_ind), 'LineStyle','none', 'Marker','*', 'MarkerSize',10, 'Color','red', 'LineWidth',2, 'HandleVisibility','off');
   if (i==1)
       legend({'$q_{d,new}$', '$q$', '$q_{new}$'}, 'interpreter','latex', 'fontsize',15);
       title('Quaternion logarithm', 'interpreter','latex', 'fontsize',17);
   end
   if (i==4), xlabel('time [$s$]', 'interpreter','latex', 'fontsize',15); end
   axis tight;
end

figure;
for i=1:3
   subplot(3,1,i); hold on;
   plot(Timed, vRotd_new_data(i,:), 'LineWidth',2, 'Color','blue', 'LineStyle','-');
   plot(Timed, vRot_data(i,:), 'LineWidth',2, 'Color',[0.85 0.33 0.1], 'LineStyle','--');
   plot(Timed, vRot_new_data(i,:), 'LineWidth',2, 'Color','magenta', 'LineStyle',':');
   plot(Timed(up_ind), vRotd_new_data(i,up_ind), 'LineStyle','none', 'Marker','*', 'MarkerSize',10, 'Color','red', 'LineWidth',2, 'HandleVisibility','off');
   if (i==1)
       legend({'$\omega_{d,new}$', '$\omega$', '$\omega_{new}$'}, 'interpreter','latex', 'fontsize',15);
       title('Rotational Velocity', 'interpreter','latex', 'fontsize',17);
   end

   if (i==4), xlabel('time [$s$]', 'interpreter','latex', 'fontsize',15); end
   axis tight;
end

figure;
for i=1:3
   subplot(3,1,i); hold on;
   plot(Timed, dvRotd_new_data(i,:), 'LineWidth',2, 'Color','blue', 'LineStyle','-');
   plot(Timed, dvRot_data(i,:), 'LineWidth',2, 'Color',[0.85 0.33 0.1], 'LineStyle','--');
   plot(Timed, dvRot_new_data(i,:), 'LineWidth',2, 'Color','magenta', 'LineStyle',':');
   plot(Timed(up_ind), dvRotd_new_data(i,up_ind), 'LineStyle','none', 'Marker','*', 'MarkerSize',10, 'Color','red', 'LineWidth',2, 'HandleVisibility','off');
   if (i==1)
       legend({'$\dot{\omega}_{d,new}$', '$\dot{\omega}$', '$\dot{\omega}_{new}$'}, 'interpreter','latex', 'fontsize',15);
       title('Rotational Acceleration', 'interpreter','latex', 'fontsize',17);
   end
   if (i==4), xlabel('time [$s$]', 'interpreter','latex', 'fontsize',15); end
   axis tight;
end
