clc;
close all;
clear;

format compact;

set_matlab_utils_path;

global Qt Rt p_o gmp_p

p_o = [0; 0; 0];
Qt = [1; 0; 0; 0];
Rt = quat2rotm(Qt');


fid = FileIO('../data/ur_train_data.bin', FileIO.in);
% fid.printHeader();

Time = fid.read('Timed');
Pd_data = fid.read('Pd_data');
dPd_data = fid.read('dPd_data');
Qd_data = fid.read('Qd_data');
% vRotd_data = fid.read('vRotd_data');

P0 = Pd_data(:,1);
Pg = Pd_data(:,end);
Q0 = Qd_data(:,1);
Qg = Qd_data(:,end);

gmp_p = GMP_nDoF(3, 30, 40, 150, 1);

setTargetPose(Pg, Qg);

n_data = size(Pd_data,2);
Pd_t_data = zeros(3,n_data);
for j=1:n_data, Pd_t_data(:,j) = w2t_pos(Pd_data(:,j)); end

setInitialPose(P0, Q0);
setTargetPose(Pg, Qg);
train_err = gmp_p.train('LS', Time, Pd_t_data);
train_error = train_err'

% simulate

ks = 1.5;
Pg = ks*(Pg-P0) + P0;

tau = Time(end) - Time(1);
x = Time/tau;
x_dot = 1/tau;
dt = Time(2) - Time(1);

Time = 0:dt:tau;
n_data = length(Time);
P_data = zeros(3, n_data);
dP_data = zeros(3, n_data);
P_t_data = zeros(3, n_data);
dP_t_data = zeros(3, n_data);


setTargetPose(Pg, Qg);
setInitialPose(P0, Q0); % do this after Pg is set!

for j=1:n_data
    %x = Time(j)/tau;
    P_data(:,j) = getRefPos(x(j));
    dP_data(:,j) = getRefVel(x(j), x_dot);
    P_t_data(:,j) = gmp_p.getYd(x(j));
    dP_t_data(:,j) = gmp_p.getYdDot(x(j), x_dot);
end

ax_legend = {{'$x$','$x_d$'}, {'$y$','$y_d$'}, {'$z$','$z_d$'}};

% fig = figure;
% fig.Position(3:4) = [1444 865];
% for i=1:3
%     subplot(2,3,i); hold on;
%     plot(Time, P_data(i,:), 'LineWidth',2, 'Color','blue');
%     plot(Time, Pd_data(i,:), 'LineWidth',2, 'Color','magenta', 'LineStyle','--');
%     legend(ax_legend{i}, 'interpreter','latex', 'fontsize',15);
%     if (i==1), ylabel('World pos [$m$]', 'interpreter','latex', 'fontsize',15); end
%         
%     subplot(2,3,i+3); hold on;
%     plot(Time, P_t_data(i,:), 'LineWidth',2, 'Color','blue');
%     plot(Time, Pd_t_data(i,:), 'LineWidth',2, 'Color','magenta', 'LineStyle','--');
%     legend(ax_legend{i}, 'interpreter','latex', 'fontsize',15);
%     xlabel('Time [$s$]', 'interpreter','latex', 'fontsize',15);
%     if (i==1), ylabel('Target pos [$m$]', 'interpreter','latex', 'fontsize',15); end    
% end

figure; hold on;
plot3(P_data(1,:), P_data(2,:), P_data(3,:), 'LineWidth',2, 'LineStyle','-', 'Color',[0.85 0.33 0.1]);
plot3(Pd_data(1,:), Pd_data(2,:), Pd_data(3,:), 'LineWidth',2, 'LineStyle','--', 'Color','blue');
scatter3(P0(1), P0(2), P0(3), 'LineWidth',4, 'SizeData',100, 'Marker','o' , 'MarkerEdgeColor','green');
scatter3(Pg(1), Pg(2), Pg(3), 'LineWidth',4, 'SizeData',100, 'Marker','x' , 'MarkerEdgeColor','red');
legend({'$p$','$p_d$', '$p_0$', '$p_g$'}, 'interpreter','latex', 'fontsize',15);

%% ===================================================================

function setTargetPose(Pg, Qg)

    global Qt Rt p_o gmp_p
    Qt = math_.quatInv(Qg);
    Rt = quat2rotm(Qt');
    p_o = Pg;
    
    gmp_p.setGoal(w2t_pos(Pg));
  
end

function setInitialPose(P0, Q0)
    
    global gmp_p
    gmp_p.setY0(w2t_pos(P0));
  
end

function p_w = getRefPos(x)

    global gmp_p
    p_w = t2w_pos(gmp_p.getYd(x));
  
end

function vel_w = getRefVel(x, x_dot)
    
    global gmp_p
	vel_w = t2w_vel(gmp_p.getYdDot(x, x_dot));
  
end

function a_w = calcAccel(pos, vel, x, x_dot, x_ddot)
    
    global gmp_p
    a_w = t2w_vel(gmp_p.calcAccel(w2t_pos(pos), w2t_vel(vel), x, x_dot, x_ddot));
  
end

function p_t = w2t_pos(pos) 

    global Rt p_o
    p_t = Rt*(pos - p_o);
    
end

function vel_t = w2t_vel(vel)
    
    global Rt
    vel_t = Rt*vel;
  
end

function p_w = t2w_pos(pos)
    
    global Rt p_o
    p_w = Rt'*pos + p_o;
  
end

function vel_w = t2w_vel(vel)
    
    global Rt
    vel_w = Rt'*vel;
  
end

