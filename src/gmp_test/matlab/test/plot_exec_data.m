clc;
% close all;
clear;

set_matlab_utils_path();

Mt_data = loadData('data/exec_data2.bin');
Mt_data2 = loadData('data/exec_data_no_perturb.bin');

%% =============================================================

k = 3; 

% [gmp_p, gmp_o] = loadModel('data/gmp_model.bin');

t_switch = Mt_data{k}.for.Time(end);
t_switch = [t_switch t_switch];
Time1 = [Mt_data{k}.for.Time Mt_data{k}.rev.Time];
phase_var = [Mt_data{k}.for.x_data Mt_data{k}.rev.x_data];
Fext = [Mt_data{k}.for.Fext Mt_data{k}.rev.Fext];
dist_data = zeros(size(Time1));
for j=1:length(dist_data), dist_data(j) = norm(Fext(:,j)); end

fig = figure('Position', [100 100 200 200]);
fig.Position(3:4) = [560 576];
% --------------------------------------------------------
% y_labels = {'x [$m$]', 'y [$m$]', 'z [$m$]'};
% for dim=1:3
%     Timed = Mt_data2{k}.for.Time;
%     Timed = Timed - Timed(1);
%     Pd_data = Mt_data2{k}.for.Pos(dim,:);
%     Timed_rev = Mt_data2{k}.rev.Time;
%     Timed_rev = Timed_rev - Timed_rev(1);
%     Pd_rev_data = Mt_data2{k}.rev.Pos(dim,:);
%     
%     Pos = [Mt_data{k}.for.Pos(dim,:) Mt_data{k}.rev.Pos(dim,:)];
% 
%     ax = subplot(4,1,dim); hold on; ax.FontSize=13;
%     plot(Time1, Pos, 'LineWidth',2.5, 'Color','blue');
%     plot(Timed+Time1(1), Pd_data, 'LineWidth',3, 'Color','magenta', 'LineStyle',':');
%     plot(Timed_rev+t_switch(1), Pd_rev_data, 'LineWidth',3, 'Color','magenta', 'LineStyle',':');
%     axis tight;
%     plot(t_switch, ax.YLim, 'LineWidth',2, 'Color',[0.85 0.33 0.1], 'LineStyle','--');
%     ylabel(y_labels{dim}, 'interpreter','latex', 'fontsize',16);
%     
%     if (dim == 2)
%        legend({'perturbed','unperturbed'}, 'interpreter','latex', 'fontsize',14, 'Position',[0.6719 0.6641 0.2727 0.0967]); 
%     end
% end
% -----------------------------------------------------------
y_labels = {'$\eta_x$ [$m$]', '$\eta_y$ [$m$]', '$\eta_z$ [$m$]'};
for dim=1:3
    Timed = Mt_data2{k}.for.Time;
    Timed = Timed - Timed(1);
    qd_data = Mt_data2{k}.for.qlog(dim,:);
    Timed_rev = Mt_data2{k}.rev.Time;
    Timed_rev = Timed_rev - Timed_rev(1);
    qd_rev_data = Mt_data2{k}.rev.qlog(dim,:);
    
    qlog = [Mt_data{k}.for.qlog(dim,:) Mt_data{k}.rev.qlog(dim,:)];

    ax = subplot(4,1,dim); hold on; ax.FontSize=13;
    plot(Time1, qlog, 'LineWidth',2.5, 'Color','blue');
    plot(Timed+Time1(1), qd_data, 'LineWidth',3, 'Color','magenta', 'LineStyle',':');
    plot(Timed_rev+t_switch(1), qd_rev_data, 'LineWidth',3, 'Color','magenta', 'LineStyle',':');
    axis tight;
    plot(t_switch, ax.YLim, 'LineWidth',2, 'Color',[0.85 0.33 0.1], 'LineStyle','--');
    ylabel(y_labels{dim}, 'interpreter','latex', 'fontsize',16);
    
    if (dim == 2)
       legend({'perturbed','unperturbed'}, 'interpreter','latex', 'fontsize',14, 'Position',[0.7237 0.6728 0.2727 0.0967]); 
    end
end
% --------------------------------------------------------
% ax = subplot(3,1,2); hold on; ax.FontSize=13;
% plot(Time1, phase_var, 'LineWidth',2.5, 'Color','blue');
% plot(Timed+Time1(1), xd_data, 'LineWidth',3, 'Color','magenta', 'LineStyle',':');
% plot(Timed_rev+t_switch(1), xd_rev_data, 'LineWidth',3, 'Color','magenta', 'LineStyle',':');
% axis tight;
% plot(t_switch, ax.YLim, 'LineWidth',2, 'Color',[0.85 0.33 0.1], 'LineStyle','--');
% ylabel('phase var', 'interpreter','latex', 'fontsize',16);
% legend({'nominal','actual'}, 'interpreter','latex', 'fontsize',17, 'Position',[0.7107 0.5433 0.2065 0.0967]);
% --------------------------------------------------------
ax = subplot(4,1,4); hold on; ax.FontSize=13;
plot(Time1, dist_data, 'LineWidth',2.5, 'Color','red');
axis tight;
plot(t_switch, ax.YLim, 'LineWidth',2, 'Color',[0.85 0.33 0.1], 'LineStyle','--');
ylabel('$d(t) = || \mathbf{F}_{ext} ||$', 'interpreter','latex', 'fontsize',16);
legend({'$d(t)$'}, 'interpreter','latex', 'fontsize',17, 'Position',[0.7790 0.2432 0.1364 0.0514]);
xlabel('time [$s$]', 'interpreter','latex', 'fontsize',16);

forward_text = annotation('textbox', 'Color',[0 0 0], 'BackgroundColor',[1 1 1], 'LineStyle','none', 'HorizontalAlignment','center',...
    'FontSize',17, 'String','Forward', 'Position',[0.2315 0.9324 0.2088 0.0521]);

reverse_text = annotation('textbox', 'Color',[1 1 1], 'BackgroundColor',[0.3137 0.3137 0.3137], 'LineStyle','none', 'HorizontalAlignment','center',...
    'FontSize',17, 'String','Reverse', 'Position',[0.6370 0.9346 0.2055 0.0587]);


return 
%% ========  Plot Cartesian Position  ==========
fig = figure;
fig.Position(3:4) = [560 420];
ax = axes();
ax.FontSize = 13;
hold on;
% -------- make legend ----------
plot3(nan,nan,nan, 'LineWidth',3, 'LineStyle','-', 'Color','green', 'DisplayName','forward');
plot3(nan,nan,nan, 'LineWidth',4, 'LineStyle',':', 'Color','red', 'DisplayName','reverse');
scatter3(nan,nan,nan, 'LineWidth',4, 'Marker','*', 'MarkerEdgeColor','magenta', 'SizeData',200, 'DisplayName','target poses');
scatter3(nan,nan,nan, 'LineWidth',4, 'Marker','o', 'MarkerEdgeColor','cyan', 'SizeData',200, 'DisplayName','initial pose');
% -------- plot ----------
for i=1:length(Mt_data)
    plot3(Mt_data{i}.for.Pos(1,:), Mt_data{i}.for.Pos(2,:), Mt_data{i}.for.Pos(3,:), 'LineWidth',3, 'LineStyle','-', 'Color','green', 'HandleVisibility','off');
    plot3(Mt_data{i}.rev.Pos(1,:), Mt_data{i}.rev.Pos(2,:), Mt_data{i}.rev.Pos(3,:), 'LineWidth',4, 'LineStyle',':', 'Color','red', 'HandleVisibility','off');
    Pg = [Mt_data{i}.for.Pos(1,end); Mt_data{i}.for.Pos(2,end); Mt_data{i}.for.Pos(3,end)];
    scatter3(Pg(1), Pg(2), Pg(3), 'LineWidth',4, 'Marker','*', 'MarkerEdgeColor','magenta', 'SizeData',200, 'HandleVisibility','off');
end
P0 = [Mt_data{1}.for.Pos(1,1); Mt_data{1}.for.Pos(2,1); Mt_data{1}.for.Pos(3,1)];
scatter3(P0(1), P0(2), P0(3), 'LineWidth',4, 'Marker','o', 'MarkerEdgeColor','cyan', 'SizeData',200, 'HandleVisibility','off');
xlabel('x [$m$]', 'interpreter','latex', 'fontsize',15);
ylabel('y [$m$]', 'interpreter','latex', 'fontsize',15);
zlabel('z [$m$]', 'interpreter','latex', 'fontsize',15);
title('Cartesian Position', 'interpreter','latex', 'fontsize',17);
legend({}, 'interpreter','latex', 'fontsize',15);
grid on;
ax.CameraPosition = [3.1267 -4.8423 3.1143];
ax.Position = [0.1300 0.1100 0.7750 0.8150];
ax.XLim = [-0.9351 0.0649];
ax.YLim = [-0.6873 0.3127];
ax.ZLim = [0.3513 0.7513];


%% ========  Plot Cartesian Orientation  ==========
fig = figure;
fig.Position(3:4) = [560 420];
ax = axes();
ax.FontSize = 13;
hold on;
% -------- make legend ----------
plot3(nan,nan,nan, 'LineWidth',3, 'LineStyle','-', 'Color','green', 'DisplayName','forward');
plot3(nan,nan,nan, 'LineWidth',4, 'LineStyle',':', 'Color','red', 'DisplayName','reverse');
scatter3(nan,nan,nan, 'LineWidth',4, 'Marker','*', 'MarkerEdgeColor','magenta', 'SizeData',200, 'DisplayName','target poses');
scatter3(nan,nan,nan, 'LineWidth',4, 'Marker','o', 'MarkerEdgeColor','cyan', 'SizeData',200, 'DisplayName','initial pose');
% -------- plot ----------
for i=1:length(Mt_data)
    plot3(Mt_data{i}.for.qlog(1,:), Mt_data{i}.for.qlog(2,:), Mt_data{i}.for.qlog(3,:), 'LineWidth',3, 'LineStyle','-', 'Color','green', 'HandleVisibility','off');
    plot3(Mt_data{i}.rev.qlog(1,:), Mt_data{i}.rev.qlog(2,:), Mt_data{i}.rev.qlog(3,:), 'LineWidth',4, 'LineStyle',':', 'Color','red', 'HandleVisibility','off');
    Pg = [Mt_data{i}.for.qlog(1,end); Mt_data{i}.for.qlog(2,end); Mt_data{i}.for.qlog(3,end)];
    scatter3(Pg(1), Pg(2), Pg(3), 'LineWidth',4, 'Marker','*', 'MarkerEdgeColor','magenta', 'SizeData',200, 'HandleVisibility','off');
end
P0 = [Mt_data{1}.for.qlog(1,1); Mt_data{1}.for.qlog(2,1); Mt_data{1}.for.qlog(3,1)];
scatter3(P0(1), P0(2), P0(3), 'LineWidth',4, 'Marker','o', 'MarkerEdgeColor','cyan', 'SizeData',200, 'HandleVisibility','off');
xlabel('x [$m$]', 'interpreter','latex', 'fontsize',15);
ylabel('y [$m$]', 'interpreter','latex', 'fontsize',15);
zlabel('z [$m$]', 'interpreter','latex', 'fontsize',15);
title('Cartesian Orientation: $\eta = \log(Q*\bar{Q}_0)$', 'interpreter','latex', 'fontsize',17);
legend({}, 'interpreter','latex', 'fontsize',15);
grid on;
ax.XLim = [-2.0000 2];
ax.YLim = [6.6613e-16 2];
ax.ZLim = [-1.5000 1.5000];
view(ax, 74.9805, 30.9473);


return
%% ==============================================================


% figure;
% ax = axes();
% hold on;
% plot(Time, P_data(1,:), 'LineWidth',2, 'Color','red');
% plot(Time, P_data(2,:), 'LineWidth',2, 'Color','green');
% plot(Time, P_data(3,:), 'LineWidth',2, 'Color','blue');
% for i=1:length(mf_ind)
%    k = mf_ind(i);
%    t1 = Time(k);
%    plot([t1 t1], ax.YLim, 'LineWidth',2, 'LineStyle','--', 'Color','magenta');
% end
% xlabel('time [$s$]', 'interpreter','latex', 'fontsize',15);
% axis tight;
% hold off;

% ------------------------------------------------------------

function Mt_data = loadData(path)

    fid = FileIO(path, FileIO.in);

    fid.printHeader();

    Time = fid.read('Time');
    x_data = fid.read('x_data');
    P_data = fid.read('P_data');
    dP_data = fid.read('dP_data');
    Q_data = fid.read('Q_data');
    vRot_data = fid.read('vRot_data');
    Fext_data = fid.read('Fext_data');
    mf_ind = fid.read('motion_finish_ind');
    
    mf_ind = mf_ind + 1;
    mf_ind = [1 mf_ind length(Time)];

    n_targets = ( length(mf_ind) - 1 ) / 2;

    Q0 = Q_data(:,1);

    Mt_data = cell(n_targets, 1);
    k = 1;
    for i=1:n_targets
        Mt_data{i} = struct('for',[], 'rev',[]);

        i1 = mf_ind(k);
        i2 = mf_ind(k+1)-1;
        Mt_data{i}.for = getSegment(Time(i1:i2), x_data(i1:i2), P_data(:, i1:i2), Q_data(:, i1:i2), Fext_data(:, i1:i2), Q0);
        k = k + 1;

        i1 = mf_ind(k);
        i2 = mf_ind(k+1)-1;
        Mt_data{i}.rev = getSegment(Time(i1:i2), x_data(i1:i2), P_data(:, i1:i2), Q_data(:, i1:i2), Fext_data(:, i1:i2), Q0);
        k = k + 1;
    end

end

% ------------------------------------------------------------

function seg = getSegment(Time, x_data, P_data, Q_data, Fext_data, Q0)
    
    qlog = zeros(3, size(Q_data,2));
    for j=1:size(qlog,2), qlog(:,j) = math_.quatLog(math_.quatDiff(Q_data(:,j),Q0)); end
    
    seg = struct('Time',Time, 'x_data',x_data, 'Pos',P_data, 'Quat',Q_data, 'qlog',qlog, 'Fext',Fext_data);

end

% ------------------------------------------------------------

function [gmp_p, gmp_o] = loadModel(path)
    
    fid = FileIO(path, FileIO.in);
    
    gmp_p = GMP_nDoF(3, 2, 1, 1);
    gmp_o = GMPo(2, 1, 1);
    
    gmp_p.readFromFile(fid, 'pos_');
    gmp_o.readFromFile(fid, 'orient_');

end

% ------------------------------------------------------------

function [Time, x_data, Y_data] = simulateGMP(gmp, y0, g, T, dt, reverse)

%% set initial values
Dim = 1;
ddy = zeros(Dim,1);
dy = zeros(Dim,1);
y = y0;
t = 0.0;
dz = zeros(Dim,1);
z = zeros(Dim,1);

t_end = T;
tau = t_end;
x = 0.0;
x_dot = 1/tau;
s = [x; x_dot];

iters = 0;
Time = [];
Y_data = [];
dY_data = [];
ddY_data = [];
x_data = [];

gmp.setY0(y0);
gmp.setGoal(g);

target = g;
dir = 1;
if (reverse)
    dir = -1;
    y = g;
    x = 1;
    target = y0;
end

%% simulate
while (true)

    %% data logging
    Time = [Time t];
    Y_data = [Y_data y];
    dY_data = [dY_data dy];  
    ddY_data = [ddY_data ddy];
    x_data = [x_data x];

    %% DMP simulation
    y_c = 0; 
    z_c = 0.0;
    gmp.update(s, y, z, y_c, z_c);
    dy = gmp.getYdot();
    dz = gmp.getZdot();

    yc_dot = 0.0;
    ddy = gmp.getYddot(yc_dot);

    %% Update phase variable
    dx = dir * 1/tau;

    %% Stopping criteria
    if (t>=1*t_end) % && norm(y-target)<5e-3 && norm(dy)<5e-3)
        break;
    end

    %% Numerical integration
    iters = iters + 1;
    t = t + dt;
    x = x + dx*dt;
    y = y + dy*dt;
    z = z + dz*dt;

    s = [x; dx];
    
end


end

% ------------------------------------------------------------