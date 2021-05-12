function plotAnimatedFrame(ax, Pos, Quat, sample_step, time_step, frame_scale)

if (nargin < 4), frame_scale = 0.1; end

prev_hold_status = ax.NextPlot;
hold(ax, 'on');

Pos = Pos(:, 1:sample_step:end);
Quat = Quat(:, 1:sample_step:end);

n_data = size(Pos, 2);

Axang = quat2axang(Quat')';
X = Pos(1,:);   Y = Pos(2,:);    Z = Pos(3,:);
U = Axang(1,:);  V = Axang(2,:);  W = Axang(3,:);  Theta = Axang(4,:);

axis_colors = {[1 0 0], [0 1 0], [0 0 1]};

% initialize frame
quiv = cell(3,1); % three quivers, for x, y and z axis of orientation frame
for i=1:3
    quiv{i} = quiver3(ax, nan,nan,nan, nan,nan,nan, frame_scale);
    quiv{i}.Color = axis_colors{i};
    quiv{i}.LineStyle = '-';
    quiv{i}.LineWidth = 2;    
    quiv{i}.AutoScale = 'on';
    quiv{i}.HandleVisibility = 'off';
end


for j=1:n_data
    
    T = makehgtform('translate',[X(j) Y(j) Z(j)]) * makehgtform('axisrotate',[U(j) V(j) W(j)],Theta(j));

    for i=1:3
        quiv{i}.XData = X(j);
        quiv{i}.YData = Y(j);
        quiv{i}.ZData = Z(j);

        quiv{i}.UData = T(1,i);
        quiv{i}.VData = T(2,i);
        quiv{i}.WData = T(3,i);
    end

    drawnow;
    pause(time_step);

%     if (capture_video)
%         frame = getframe(ax);
%         videoWriter.writeVideo(frame);
%     end

end

for i=1:3, delete(quiv{i}); end

ax.NextPlot = prev_hold_status;

end