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

tf = 2;
dt = 0.001;
n = round(tf/dt);

Time = zeros(1,n);

t = 0;

y = y0;
y_dot = 0;
y_ddot = 0;
z = y_dot;
z_dot = 0;
s = 0;
s_dot = 0;
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
fp_data = zeros(1,n);

a1 = 1e-8;
a2 = 0;

% potential_fun = @simple_potential;
potential_fun = @ppc_potential;

yL_up = @(t) 0.8*ones(size(t)) + 0.2*sin(2*pi*t/1);
yL_low = @(t) 0.2*ones(size(t)) + 0.1*sin(2*pi*t/1);
yL_dot_up = @(t) 1.5*ones(size(t)) + 0.3*sin(2*pi*t/1);
yL_dot_low = @(t) -1.5*ones(size(t)) + 0.3*sin(2*pi*t/1);

K = 300;
D = 30;

Kp = 0;

a_pd = 1e-8;
a_p = a_pd;

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
    yd_ddot = (fs_dot(t+dt)-fs_dot(t))/dt;
    
%     a_p = 1e-5;
%     fp = a_p*( 1/(y - yL_up(t))^3 + 1/(y - yL_low(t))^3);
    fp = simple_potential(y, yL_low(t), yL_up(t), a_p);
%     fp = sign(fp)*min([fp 10]);
%     fp = ppc_potential(y, yL_low(t), yL_up(t), 10);
%     fp_dot = -3*a_p*( 1/(y - yL_up(t))^4 + 1/(y - yL_low(t))^4 )*y_dot;
%     fv = 40*( 1/(y_dot - yL_dot_up(t))^3 + 1/(y_dot - yL_dot_low(t))^3);
    fv = simple_potential(y_dot, yL_dot_low(t), yL_dot_up(t), 40);
    
%     y_ddot = K*(yd-y) + D*(yd_dot + fp - y_dot) + fv; % no use of fp_dot
%     y_ddot = K*(yd-y) + D*(yd_dot + fp - y_dot) + fp_dot + fv; % use fp_dot
    
    flag = true;
    y_low = yL_low(t);
    y_up = yL_up(t);
%     while (flag)
        z_dot = K*(yd-y) + D*(yd_dot - z) + 0*fv; % use fp_dot indirectly
        y_dot = z + fp;
    
%         break;
        
        y_hat = y + y_dot*dt;
        z_hat = z + z_dot*dt;
        
        
%         if ( y_hat > y_low && y_hat < y_up )
%             flag = false;
%         else
%             yi = y_low + 1e-3;
%             if (norm(y-y_up) < norm(y-yi)), yi = y_up-1e-3; end
% %             y_low
% %             y_up
% %             y_hat
% %             y_dot
% %             z
% %             y + (z + fp)*dt
%             fp = (yi-y)/dt - z;
%             a_p = fp / simple_potential(y, yL_low(t), yL_up(t), 1);
%             t
% %             y_hat = y + (z + fp)*dt
% %             yi
% %             pause
%         end
%     end
    
%     z_dot = K*(yd-y) + D*(yd_dot - z) + fp;
%     y_dot = z;

    fp_data(j) = fp;
    
    y2_ddot = K*(yd-y2) + D*(yd_dot - y2_dot);
    
    t = t + dt;
    y = y + y_dot*dt;
    z = z + z_dot*dt;
    %s = s + s_dot*dt;
    %y_dot = y_dot + y_ddot*dt;
    
    y2 = y2 + y2_dot*dt;
    y2_dot = y2_dot + y2_ddot*dt;
    
end


% figure;
% subplot(2,1,1); hold on;
% plot(Time, y_data, 'LineWidth',2);
% plot(Time, y2_data, 'LineWidth',2, 'Color',[0 0.6 0 0.5], 'LineStyle',':');
% plot(Time, yL_up(Time), 'LineWidth',2, 'Color','red', 'LineStyle','--');
% plot(Time, yL_low(Time), 'LineWidth',2, 'Color','magenta', 'LineStyle','--');
% subplot(2,1,2); hold on;
% plot(Time, fp_data, 'LineWidth',2);

figure;
subplot(2,1,1); hold on;
plot(Time, y_data, 'LineWidth',2);
plot(Time, y2_data, 'LineWidth',2, 'Color',[0 0.6 0 0.5], 'LineStyle',':');
plot(Time, yL_up(Time), 'LineWidth',2, 'Color','red', 'LineStyle','--');
plot(Time, yL_low(Time), 'LineWidth',2, 'Color','magenta', 'LineStyle','--');
subplot(2,1,2); hold on;
plot(Time, y_dot_data, 'LineWidth',2);
plot(Time, y2_dot_data, 'LineWidth',2, 'Color',[0 0.6 0 0.5], 'LineStyle',':');
plot(Time, yL_dot_up(Time), 'LineWidth',2, 'Color','red', 'LineStyle','--');
plot(Time, yL_dot_low(Time), 'LineWidth',2, 'Color','magenta', 'LineStyle','--');


function fp = simple_potential(y, y_low, y_up, gain)

    fp = gain * ( 1/(y - y_up)^3 + 1/(y - y_low)^3);
%     fp = gain * 1 /( (y - y_up)*(y - y_low) );

end

function fp = ppc_potential(y, y_low, y_up, gain)

    fp = gain * log10((y-y_low)/(y_up-y)) * (y_up - y_low) / ( (y - y_low)*(y_up - y) );

end

