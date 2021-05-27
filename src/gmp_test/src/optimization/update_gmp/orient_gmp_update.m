% function orient_gmp_test()

clc;
% close all;
clear;

addpath('../../../../matlab/lib/gmp_lib/');
import_gmp_lib();

gmp_o = GMPo();
gmp_.read(gmp_o, 'data/gmp_orient.bin','');

Ts = 0.01;
Timed = 0:Ts:10;
x = Timed/Timed(end);
n_data = length(x);

%% Construct new orientation trajectory

% set demo and new init/target orientation
Q0d = gmp_o.getQd(0);
Qgd = gmp_o.getQd(1);

Q0 = Q0d;
Qg = gmp_.quatProd( gmp_.quatExp([0.5, 0.2, -0.4]), Qgd);

% generate nominal trajectory for the new init/target orientation
gmp_o.setQ0(Q0);
gmp_o.setQg(Qg);
for j=1:n_data, Qd_data(:,j) = gmp_o.getQd(x(j)); end


qd_data = zeros(3, n_data);

Qd_new = zeros(4, n_data);
qd_new_data = zeros(3, n_data);
vRotd_new = zeros(3, n_data);
dvRotd_new = zeros(3, n_data);

for j=1:n_data
    qd_data(:,j) = GMPo.quat2q(Qd_data(:,j), Q0);
    qd_new_data(:,j) = qd_data(:,j) + 0.3*exp(-((0.4-x(j))/0.1)^2) - 0.2*exp(-((0.7-x(j))/0.06)^2);
    Qd_new(:,j) = GMPo.q2quat(qd_new_data(:,j), Q0);
end

for j=1:n_data-1
    vRotd_new(:,j) = math_.quatLog( math_.quatDiff(Qd_new(:,j+1),Qd_new(:,j)) )/Ts;
end
for i=1:3, dvRotd_new(i,:)=[diff(vRotd_new(i,:)) 0]/Ts; end


%% Update GMP

tau = Timed(end);
x_dot = 1/tau;
x_ddot = 0;

gmp_o_new = gmp_o.deepCopy();
gmp_up = GMPo_Update(gmp_o_new);
gmp_up.enableSigmawUpdate(true);
gmp_up.setMsrNoiseVar(1e-4);

up_ind = [];
t_span = 0;
for j=1:n_data
    t_span = t_span + 1000*Ts;
    if ( norm(qd_new_data(:,j)-qd_data(:,j) ) > 0.05  && t_span>300) % update every 300 ms
        up_ind = [up_ind j];
        t_span = 0;
    end
end
x_up = x(up_ind);


for j=1:length(up_ind)
    k = up_ind(j);
    gmp_up.updatePos(x_up(j), qd_new_data(:,k));
%     gmp_up.updateQuat(x_up(j), Qd_new(:,k));
%     gmp_up.updateRotVel(x_up(j), x_dot, vRotd_new(:,k), Qd_new(:,k) );
%     gmp_up.updateRotAccel(x_up(j), x_dot, x_ddot, dvRotd_new(:,k), vRotd_new(:,k), Qd_new(:,k) );
end

q_data = zeros(3, n_data);
Q_data = zeros(4, n_data);
vRot_data = zeros(3, n_data);
dvRot_data = zeros(3, n_data);

q_new_data = zeros(3, n_data);
Q_new_data = zeros(4, n_data);
vRot_new_data = zeros(3, n_data);
dvRot_new_data = zeros(3, n_data);

for j=1:n_data  
    q_data(:,j) = gmp_o.getYd(x(j));
    Q_data(:,j) = gmp_o.getQd(x(j));
    vRot_data(:,j) = gmp_o.getVd(x(j), x_dot);
    dvRot_data(:,j) = gmp_o.getVdDot(x(j), x_dot, x_ddot);

    q_new_data(:,j) = gmp_o_new.getYd(x(j));
    Q_new_data(:,j) = gmp_o_new.getQd(x(j));
    vRot_new_data(:,j) = gmp_o_new.getVd(x(j), x_dot);
    dvRot_new_data(:,j) = gmp_o_new.getVdDot(x(j), x_dot, x_ddot);
end


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
   plot(Timed, vRotd_new(i,:), 'LineWidth',2, 'Color','blue', 'LineStyle','-');
   plot(Timed, vRot_data(i,:), 'LineWidth',2, 'Color',[0.85 0.33 0.1], 'LineStyle','--');
   plot(Timed, vRot_new_data(i,:), 'LineWidth',2, 'Color','magenta', 'LineStyle',':');
   plot(Timed(up_ind), vRotd_new(i,up_ind), 'LineStyle','none', 'Marker','*', 'MarkerSize',10, 'Color','red', 'LineWidth',2, 'HandleVisibility','off');
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
   plot(Timed, dvRotd_new(i,:), 'LineWidth',2, 'Color','blue', 'LineStyle','-');
   plot(Timed, dvRot_data(i,:), 'LineWidth',2, 'Color',[0.85 0.33 0.1], 'LineStyle','--');
   plot(Timed, dvRot_new_data(i,:), 'LineWidth',2, 'Color','magenta', 'LineStyle',':');
   plot(Timed(up_ind), dvRotd_new(i,up_ind), 'LineStyle','none', 'Marker','*', 'MarkerSize',10, 'Color','red', 'LineWidth',2, 'HandleVisibility','off');
   if (i==1)
       legend({'$\dot{\omega}_{d,new}$', '$\dot{\omega}$', '$\dot{\omega}_{new}$'}, 'interpreter','latex', 'fontsize',15);
       title('Rotational Acceleration', 'interpreter','latex', 'fontsize',17);
   end
   if (i==4), xlabel('time [$s$]', 'interpreter','latex', 'fontsize',15); end
   axis tight;
end
