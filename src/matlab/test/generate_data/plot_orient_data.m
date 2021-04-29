clc;
close all;
clear;

addpath('utils/');

load('data/orient_data.mat','Data');

figure;
for i=1:4
   subplot(3,12, [1 2 3] + (i-1)*3 );
   plot(Data.Time, Data.Quat(i,:), 'LineWidth',2, 'LineStyle','-', 'Color','blue');
   if (i==1), ylabel('$Q$', 'interpreter','latex', 'fontsize',15);  end
end
for i=1:3
   subplot(3,12, 12 + [1 2 3 4] + (i-1)*4 );
   plot(Data.Time, Data.rotVel(i,:), 'LineWidth',2, 'LineStyle','-', 'Color','green');
   if (i==1), ylabel('$\omega$', 'interpreter','latex', 'fontsize',15);  end
end
for i=1:3
   subplot(3,12, 24 + [1 2 3 4] + (i-1)*4 );
   plot(Data.Time, Data.rotAccel(i,:), 'LineWidth',2, 'LineStyle','-', 'Color','red');
   if (i==1), ylabel('$\dot{\omega}$', 'interpreter','latex', 'fontsize',15);  end
end


