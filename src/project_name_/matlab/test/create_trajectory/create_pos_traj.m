clc;
close all;
clear;

p0 = 0.45*[0.8 0.6 0.9]';

a = (2.4/norm(p0))^1.3*[2; 2; 0.6];

Time = [];
P_data = [];
dP_data = [];

dt = 0.002;
t = 0;
p = p0;
p_dot = zeros(3,1);

while (norm(p) > 1e-3)
    
    Time = [Time t];
    P_data = [P_data p];
    dP_data = [dP_data p_dot];
    
    a1 = a*(1-exp(-0.1*t))*exp(0.05*t);
    p_dot = -a1.*p;
    
    t = t + dt;
    p = p + p_dot*dt;
    
end

load('data.mat', 'Data');
Data.Time = Time;
Data.Pos = P_data;
save('data.mat', 'Data');

figure;
for i=1:3
    subplot(3,1,i);
    plot(Time, dP_data(i,:), 'LineWidth',2);
    axis tight;
    if (i==3), xlabel('time [$s$]', 'interpreter','latex', 'fontsize',14); end
end

figure; hold on;
plot3(P_data(1,:), P_data(2,:), P_data(3,:), 'LineWidth',2, 'Color',[0.85 0.33 0.1]);
scatter3(p0(1), p0(2), p0(3), 'LineWidth',4, 'SizeData',150, 'Marker','o', 'MarkerEdgeColor','green');
scatter3(0, 0, 0, 'LineWidth',4, 'SizeData',150, 'Marker','x', 'MarkerEdgeColor','red');
xlabel('$x$ [$m$]', 'interpreter','latex', 'fontsize',15);
ylabel('$y$ [$m$]', 'interpreter','latex', 'fontsize',15);
zlabel('$z$ [$m$]', 'interpreter','latex', 'fontsize',15);
view(149.6, 22.8);
axis tight;

norm(p0)
t_end = Time(end)
