clear;
clc;
close all;


set_matlab_utils_path();

%% Load training data
% load('data/train_data.mat', 'Data');
load('data/train_data1.mat', 'Data');

y0_d = 0;
yf_d = 1;
Timed = 0:0.002:5;
[Yd_data, dYd_data, ddYd_data] = get5thOrderPol(y0_d, yf_d, Timed);

Ts = Timed(2)-Timed(1);

%% initialize DMP
N_kernels = 30;
K = 5000;
D = 2*sqrt(K);
a_z = D;
b_z = K/D;
can_clock_ptr = CanonicalClock();
shape_attr_gat_ptr = SigmoidGatingFunction(1.0, 0.5);
% shape_attr_gat_ptr = LinGatingFunction(1.0, 0.05);
% shape_attr_gat_ptr = LinGatingFunction(1.0, 0.01);
dmp = DMP(N_kernels, a_z, b_z, can_clock_ptr, shape_attr_gat_ptr);
offline_train_mse = dmp.train('LS', Timed, Yd_data, dYd_data, ddYd_data)


%% simulate DMP

rev = true; % reverse DMP flag

y0 = Yd_data(:, 1);
g = Yd_data(:, end);

dt = Ts;

x = 0.0;
dx = 0.0;
ddy = 0;
dy = 0;
y = y0;
t = 0.0;

t_end = Timed(end);
tau = t_end;
can_clock_ptr.setTau(tau);

iters = 0;
Time = [];
Y_data = [];
dY_data = [];
ddY_data = [];
x_data = [];

dmp.setY0(y0);

dir = 1;
if (rev)
    dir = -1;
    y = g;
    x = 1;
end

%% simulate
while (true)

    %% Stopping criteria
    if (t>=1.0*t_end) % && norm(y-g)<5e-3 && norm(dy)<5e-3)
        break;
    end
    
    %% data logging
    Time = [Time t];
    Y_data = [Y_data y];
    dY_data = [dY_data dy];  
    ddY_data = [ddY_data ddy];
    % x_data = [x_data x];

    %% DMP simulation
    ddy = dmp.calcYddot(x, y, dy, g);

    %% Update phase variable
    dx = dir * can_clock_ptr.getPhaseDot(x);

    %% Numerical integration
    iters = iters + 1;
    t = t + dt;
    x = x + dx*dt;
    if (x <= 0), x = 0; end
    dy = dy + ddy*dt;
    y = y + dy*dt;

end

Yd_rev_data = fliplr(Yd_data);

rmse = sqrt(mean((Y_data - Yd_rev_data).^2))
mae = mean(abs(Y_data - Yd_rev_data))

%% Plot results
figure;
ax = axes();
hold on
plot(Timed, Yd_rev_data, 'LineWidth',3, 'LineStyle','-', 'Color','cyan');
plot(Time, Y_data, 'LineWidth',4, 'LineStyle',':', 'Color','magenta');
scatter(Time(end), y0, 'LineWidth',4, 'Marker','o', 'MarkerEdgeColor',[0 1 0], 'SizeData',200);
scatter(Time(1), g, 'LineWidth',4, 'Marker','x', 'MarkerEdgeColor',[1 0 0], 'SizeData',200);
xlabel('time [$s$]', 'interpreter','latex', 'fontsize',15);
ylabel('pos [$m$]', 'interpreter','latex', 'fontsize',15);
legend({'desired','rev DMP'}, 'interpreter','latex', 'fontsize',15, 'Position',[0.6854 0.7258 0.2236 0.1250]);
axis tight;
% ax.XLim = ax.XLim + [-0.05 0.05];
% ax.YLim = ax.YLim + [-0.05 0.05];


%% ===================================================================



