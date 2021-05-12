clc;
close all;
clear;

% s = tf('s');
% H = (s+14)^3
% return

fp = @(t) sin(2*pi*t/2);
fp_dot = @(t) (2*pi/2)*cos(2*pi*t/2);

p = 0;
p_dot = 0;

t = 0;
tf = 6;
dt = 0.002;

p_lim = @(t) 0.6*ones(size(t)) + 0.25*sin(2*pi*t/1.2);
dp_lim = @(t) 1.2*ones(size(t)) + 0.3*sin(2*pi*t/1.5);

y = p; % p
z = 0; % p_dot
w = 0; % p_ddot

Time = [];
y_data = [];
z_data = [];

while (t <= tf)
    
    % logging
    Time = [Time t];
    y_data = [y_data y];
    z_data = [z_data z];
   
    % dynamic equations
    y_dot = z + 1e-3*( 1/(y - p_lim(t))^3 + 1/(y + p_lim(t))^3);
    z_dot = w + 5e-1*( 1/(z - dp_lim(t))^3 + 1/(z + dp_lim(t))^3);
    w_dot = 2744*(fp(t) - y) + 588*(fp_dot(t) - z) - 42*w;
     
    % integration
    t = t + dt;
    y = y + y_dot*dt;
    z = z + z_dot*dt;
    w = w + w_dot*dt;
    
end

fp_data = fp(Time);
fp_dot_data = fp_dot(Time);

figure('Position',[488.2 41.8 560 740.8]);
subplot(2,1,1); hold on;
plot(Time, y_data, 'LineWidth',2, 'LineStyle','-', 'Color','blue');
plot(Time, fp_data, 'LineWidth',2, 'LineStyle',':', 'Color','magenta');
plot(Time, p_lim(Time), 'LineWidth',2, 'LineStyle','--', 'Color','red');
plot(Time, -p_lim(Time), 'LineWidth',2, 'LineStyle','--', 'Color','red');
axis tight;
subplot(2,1,2); hold on;
plot(Time, z_data, 'LineWidth',2, 'LineStyle','-', 'Color','blue');
plot(Time, fp_dot_data, 'LineWidth',2, 'LineStyle',':', 'Color','magenta');
plot(Time, dp_lim(Time), 'LineWidth',2, 'LineStyle','--', 'Color','red');
plot(Time, -dp_lim(Time), 'LineWidth',2, 'LineStyle','--', 'Color','red');
axis tight;
