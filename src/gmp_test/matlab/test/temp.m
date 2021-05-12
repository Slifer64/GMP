clear;
clc;
% close all;


set_matlab_utils_path();

%% Load training data
% load('data/train_data.mat', 'Data');
load('data/train_data1.mat', 'Data');

ind = [1 3];

Timed = Data.Time;
Pd_data = Data.Pos(ind,:);
dPd_data = Data.Vel(ind,:);
ddPd_data = Data.Accel(ind,:);

Ts = Timed(2)-Timed(1);

train_method = 'LS';
N_kernels = 30;
%% initialize DMP
K = 5000;
D = 2*sqrt(K);
a_z = D;
b_z = K/D;
can_clock_ptr = CanonicalClock();
shape_attr_gat_ptr = SigmoidGatingFunction(1.0, 0.5);
% shape_attr_gat_ptr = LinGatingFunction(1.0, 0.05);
% shape_attr_gat_ptr = LinGatingFunction(1.0, 0.01);
dmp = cell(2,1);
for i=1:length(dmp)
    dmp{i} = DMP(N_kernels, a_z, b_z, can_clock_ptr, shape_attr_gat_ptr);
    offline_train_mse = dmp{i}.train(train_method, Timed, Pd_data(i,:), dPd_data(i,:), ddPd_data(i,:));
end

gmp = GMP_nDoF(size(Pd_data,1), N_kernels, D, K, 1.0);
offline_train_mse = gmp.train(train_method, Timed, Pd_data);


reverse = true;

ks = 1;
kt = 1;

y0 = Pd_data(:, 1);
Pgd = Pd_data(:, end);
g = y0 + ks*(Pgd - y0);
T = Timed(end) / kt;
dt = Ts;


Dim = length(dmp);
x = 0.0;
dx = 0.0;
ddy = zeros(Dim,1);
dy = zeros(Dim,1);
y = y0;
t = 0.0;

t_end = T;
tau = T;
can_clock_ptr.setTau(tau);

iters = 0;
Time = [];
Y_data = [];
dY_data = [];
ddY_data = [];
x_data = [];

for i=1:Dim, dmp{i}.setY0(y0(i)); end

dir = 1;
if (reverse)
    dir = -1;
    y = g;
    x = 1;
end

%% simulate
while (true)

    %% data logging
    Time = [Time t];
    Y_data = [Y_data y];
    dY_data = [dY_data dy];  
    ddY_data = [ddY_data ddy];
    % x_data = [x_data x];
    
    yd_dot = gmp.getYdDot(x, 1/tau);

    %% DMP simulation
    for i=1:Dim
        ddy(i) = dmp{i}.calcYddot(x, y(i), dy(i), g(i)) - 0*a_z*yd_dot(i)/tau;
    end

    %% Update phase variable
    dx = dir * can_clock_ptr.getPhaseDot(x);

    %% Stopping criteria
    if (t>=1.5*t_end) % && norm(y-g)<5e-3 && norm(dy)<5e-3)
        break;
    end

    %% Numerical integration
    iters = iters + 1;
    t = t + dt;
    x = x + dx*dt;
    if (x <= 0), x = 0; end
    dy = dy + ddy*dt;
    y = y + dy*dt;

end



%% Plot results
figure;
ax = axes();
hold on
plot(Pd_data(1,:), Pd_data(2,:), 'LineWidth',3 ,  'LineStyle','-', 'Color','cyan');
plot(Y_data(1,:), Y_data(2,:), 'LineWidth',4, 'LineStyle',':', 'Color','magenta');
scatter(y0(1), y0(2), 'LineWidth',4, 'Marker','o', 'MarkerEdgeColor',[0 1 0], 'SizeData',200);
scatter(g(1), g(2), 'LineWidth',4, 'Marker','x', 'MarkerEdgeColor',[1 0 0], 'SizeData',200);
xlabel('x-axis [$m$]', 'interpreter','latex', 'fontsize',15);
ylabel('y-axis [$m$]', 'interpreter','latex', 'fontsize',15);
% legend({'forward','reverse','initial','target'}, 'interpreter','latex', 'fontsize',15, 'Position',[0.6840 0.2293 0.1891 0.2256]);
axis tight;
ax.XLim = ax.XLim + [-0.05 0.05];
ax.YLim = ax.YLim + [-0.05 0.05];
hold off;


%% ===================================================================



