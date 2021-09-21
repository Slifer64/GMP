clc;
close all;
clear;

x(:,1:2) = repmat([0 0 1 1]', 1, 2);
y(:,1:2) = repmat([0 1 1 0]', 1, 2);
z(:,1:2) = [[0 0 0 0]' [1 1 1 1]'];

figure;
patch(x,y,z, 'r');
view(10, 10);

return



a = -pi : pi/2 : pi;                                % Define Corners
ph = pi/4;                                          % Define Angular Orientation (‘Phase’)
x = [cos(a+ph); cos(a+ph)]/cos(ph);
y = [sin(a+ph); sin(a+ph)]/sin(ph);
z = [-ones(size(a)); ones(size(a))];
figure
surf(x, y, z, 'FaceColor','g')                      % Plot Cube
hold on
patch(x', y', z', 'r')                              % Make Cube Appear Solid
hold off
axis([ -1  1    -1  1    -1  1]*1.5)
grid on