% clc;
% close all;
% clear;

set_matlab_utils_path();

y0 = 0;
yf = 1;

dt = 0.005;
tf = 2;
Timed = 0:dt:tf;

[yd, yd_dot, yd_ddot] = get5thOrderPol(y0, yf, Timed);


figure;
plot(Timed, yd, 'LineWidth',2);
xlabel('time [$s$]', 'interpreter','latex', 'fontsize',15);


Data = struct('Time',Timed, 'Pos',yd, 'Vel',yd_dot, 'Accel',yd_ddot);
save('data/fifthOrd_train_data.mat', 'Data');



