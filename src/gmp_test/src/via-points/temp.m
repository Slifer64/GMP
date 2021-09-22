clc;
close all;
clear;


figure;
ax = axes('NextPlot','add');

% plotCircle(ax, 1, [1 1 0], [1 0 0 pi/4]);

% [X,Y,Z] = meshSphere(2, [2 1 1]);

% [X,Y,Z] = meshCylinder(0.5, 2, [0 0 0], [1 0 0 pi/4]);

[X,Y,Z] = meshBox(4, 1, 2, [0 0 0]);

% [X,Y,Z] = meshCircle(1, [1 1 1], [1 0 0 pi/4]);

s = surf(X,Y,Z, 'EdgeColor','none', 'FaceColor',[0.8 0 0], 'FaceAlpha',0.5, 'Parent',ax);

xlabel('X', 'fontsize',15);
ylabel('Y', 'fontsize',15);
zlabel('Z', 'fontsize',15);

axis([repmat([-5 5],1,3)]);

view(-41,23);

% for i=1:70
%     pause(0.02);
%     translateSurface(s, [-0.05 -0.05 -0.01]);
%     drawnow;
% end

return;


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


function [X,Y,Z] = meshCircle(radius, center, axis_angle)
    
    if (nargin < 2), center = [0 0 0]; end
    if (nargin < 3), axis_angle = [0 0 1 0]; end
    
    [r_i, theta_i] = meshgrid(0:radius/20:radius, 0:2*pi/20:2*pi);
    X = r_i.*cos(theta_i); % + center(1);
    Y = r_i.*sin(theta_i); % + center(2); 
    Z = zeros(size(X)); % + center(3); 
    
    [m,n] = size(X);
    
    Rotm = axang2rotm(axis_angle);
    temp = Rotm*[X(:) Y(:) Z(:)]';
    X = reshape(temp(1,:), m, n) + center(1);
    Y = reshape(temp(2,:), m, n) + center(2);
    Z = reshape(temp(3,:), m, n) + center(3);
    % rotate(s, axis_angle(1:3), axis_angle(4));
    
end

function [X,Y,Z] = meshSphere(radius, center)

    if (nargin < 2), center = [0 0 0]; end

    [theta_i, phi_j] = meshgrid(0:2*pi/20:2*pi, 0:pi/20:pi);
    X = radius*cos(theta_i).*sin(phi_j) + center(1);
    Y = radius*sin(theta_i).*sin(phi_j) + center(2);
    Z = radius*cos(phi_j) + center(3);

end

function [X,Y,Z] = meshCylinder(radius, height, center, axis_angle)

    if (nargin < 3), center = [0 0 0]; end
    if (nargin < 4), axis_angle = [0 0 1 0]; end
    
    [theta, h] = meshgrid(0:2*pi/20:2*pi, [0 height]);
    X = radius*cos(theta);
    Y = radius*sin(theta);
    Z = h;
    
    [m,n] = size(X);
    
    Rotm = axang2rotm(axis_angle);
    temp = Rotm*[X(:) Y(:) Z(:)]';
    X = reshape(temp(1,:), m, n) + center(1);
    Y = reshape(temp(2,:), m, n) + center(2);
    Z = reshape(temp(3,:), m, n) + center(3);

end

function [X,Y,Z] = meshBox(x_len, y_len, z_len, center)
    
    X = (x_len/2)*repmat([-1 -1 1  1 -1],2,1);
    Y = (y_len/2)*repmat([-1  1 1 -1 -1],2,1);
    Z = (z_len/2)*[-ones(1,5); ones(1,5)];
    
end

function s = translateSurface(s, translation)

    s.XData = s.XData + translation(1);
    s.YData = s.YData + translation(2);
    s.ZData = s.ZData + translation(3);
    
end