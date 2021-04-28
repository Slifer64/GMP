clc;
close all;
clear;

g = 0.5;
y0 = 0.5;

fs = @(x) 0.5 + 0.7*sin(2*pi*x/1);
fs_dot = @(x) 0.7*(2*pi/1)*cos(2*pi*x/1);

% Time = 0:dt:5;
% f = fs(Time);
% figure;
% plot(Time, f, 'LineWidth',2);

n = 1000;
Time = zeros(1,n);

dt = 0.002;
t = 0;

y = y0;
y_dot = 0;
y_ddot = 0;
z = y_dot;
y_data = zeros(1,n);
y_dot_data = zeros(1,n);
y_ddot_data = zeros(1,n);

y2 = y0;
y2_dot = 0;
y2_ddot = 0;
z2 = y2_dot;
y2_data = zeros(1,n);
y2_dot_data = zeros(1,n);
y2_ddot_data = zeros(1,n);

a1 = 1e-5;
a2 = 0;

yL_up = 0.8;
yL_low = 0.2;
yL_dot_up = 1.5;
yL_dot_low = -1.5;

K = 300;
D = 30;

Kp = 0;

for j=1:n
    
    Time(j) = t;
    y_data(j) = y;
    y_dot_data(j) = y_dot;
    y_ddot_data(j) = y_ddot;
    
    y2_data(j) = y2;
    y2_dot_data(j) = y2_dot;
    y2_ddot_data(j) = y2_ddot;
    
    yd = fs(t);
    yd_dot = fs_dot(t);

    a = 1e-3;
    
%     z_dot = K*(yd-y) + D*(yd_dot - z); % + a2/(y_dot - yL_dot_up)^3;% + a2/(y_dot - yL_dot_low)^3;
%     y_dot = z + a/(y - yL_up)^3 + a/(y - yL_low)^3 + 0*Kp*(yd-y);
    
    
    f = a*( 1/(y - yL_up)^3 + 1/(y - yL_low)^3);
    f_dot = -3*a*( 1/(y - yL_up)^4 + 1/(y - yL_low)^4 )*y_dot;
    
    y_ddot = K*(yd-y) + D*(yd_dot + f - y_dot) + f_dot;
    
    
%     y_ddot_ = z_dot + f_dot;
%     y_ddot - y_ddot
%     pause
    
    y2_ddot = K*(yd-y2) + D*(yd_dot - y2_dot);
    
    t = t + dt;
    y = y + y_dot*dt;
%     z = z + z_dot*dt;
    y_dot = y_dot + y_ddot*dt;
    
    y2 = y2 + y2_dot*dt;
    y2_dot = y2_dot + y2_ddot*dt;
    
end


figure;
subplot(2,1,1); hold on;
plot(Time, y_data, 'LineWidth',2);
plot(Time, y2_data, 'LineWidth',2, 'Color',[0 0.6 0 0.5], 'LineStyle',':');
plot([Time(1) Time(end)], [yL_up yL_up], 'LineWidth',2, 'Color','red', 'LineStyle','--');
plot([Time(1) Time(end)], [yL_low yL_low], 'LineWidth',2, 'Color','magenta', 'LineStyle','--');
subplot(2,1,2); hold on;
plot(Time, y_dot_data, 'LineWidth',2);
plot(Time, y2_dot_data, 'LineWidth',2, 'Color',[0 0.6 0 0.5], 'LineStyle',':');
plot([Time(1) Time(end)], [yL_dot_up yL_dot_up], 'LineWidth',2, 'Color','red', 'LineStyle','--');
plot([Time(1) Time(end)], [yL_dot_low yL_dot_low], 'LineWidth',2, 'Color','magenta', 'LineStyle','--');



