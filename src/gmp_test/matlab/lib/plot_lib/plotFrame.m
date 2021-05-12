function plotFrame(ax, pos, quat, frame_scale)

if (nargin < 4), frame_scale = 0.1; end

prev_hold_status = ax.NextPlot;
hold(ax, 'on');

x = pos(1);
y = pos(2);
z = pos(3);

axang = quat2axang(quat')';
u = axang(1);
v = axang(2);
w = axang(3);
theta = axang(4);
T = makehgtform('translate',[x y z]) * makehgtform('axisrotate',[u v w],theta);

axis_colors = {[1 0 0], [0 1 0], [0 0 1]};

quiv = cell(3,1); % three quivers, for x, y and z axis of orientation frame
for j=1:3
    quiv{j} = quiver3(ax, x,y,z, T(1,j),T(2,j),T(3,j), frame_scale);
    quiv{j}.Color = axis_colors{j};
    quiv{j}.LineStyle = '-';
    quiv{j}.LineWidth = 2;    
    quiv{j}.AutoScale = 'on';
    quiv{j}.HandleVisibility = 'off';
end

ax.NextPlot = prev_hold_status;

end