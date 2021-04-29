clc;
close all;
clear;

addpath('utils/');

load('data/pos_data.mat','Data');

figure;
k=[1 4 7];
ax = [];
for i=1:3
    ax = [ax subplot(3,3,k(1))];
    hold on;
    plot(Data.Time, Data.Pos(i,:), 'LineWidth',2.0 , 'Color','blue');
    if (i==1), ylabel('pos', 'interpreter','latex', 'fontsize',15); end
    axis tight;

    ax = [ax subplot(3,3,k(2))];
    plot(Data.Time, Data.Vel(i,:), 'LineWidth',2.0, 'Color','green');
    if (i==1), ylabel('vel', 'interpreter','latex', 'fontsize',15); end
    axis tight;

    ax = [ax subplot(3,3,k(3))];
    plot(Data.Time, Data.Accel(i,:), 'LineWidth',2.0, 'Color','red');
    if (i==1), ylabel('accel', 'interpreter','latex', 'fontsize',15); end
    axis tight;
    
    k=k+1;
end
title_ = {'x', 'y', 'z'};
k=[1 4 7];
for i=1:3, title(title_{i}, 'interpreter','latex', 'fontsize',18, 'Parent',ax(k(i))); end

figure; hold on;
plot3(Data.Pos(1,:), Data.Pos(2,:), Data.Pos(3,:), 'LineWidth',2, 'LineStyle','-', 'Color','blue');
scatter3(Data.Pos(1,1), Data.Pos(2,1), Data.Pos(3,1), 'Marker','o', 'SizeData',100, 'LineWidth',4, 'MarkerEdgeColor','green');
scatter3(Data.Pos(1,end), Data.Pos(2,end), Data.Pos(3,end), 'Marker','*', 'SizeData',100, 'LineWidth',4, 'MarkerEdgeColor','red');
title('3D path', 'interpreter','latex', 'fontsize',18);
xlabel('x', 'interpreter','latex', 'fontsize',17);
ylabel('y', 'interpreter','latex', 'fontsize',17);
zlabel('z', 'interpreter','latex', 'fontsize',17);
