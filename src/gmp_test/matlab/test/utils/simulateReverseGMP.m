function [Time, Y_data, dY_data, ddY_data, F_dist, x_data] = simulateReverseGMP(gmp, y0, g, T, dt, disturbance_fun)
%% Simulates a dmp
% @param[in] gmp: Dim x 1 cell array, where each cell is a 1D GMP.
% @param[in] y0: Dim x 1 vector with the initial position..
% @param[in] g: Dim x 1 vector with the goal/target position.
% @param[in] T: Movement total time duration.
% @param[in] dt: Simulation timestep.
% @param[out] Time: 1 x N rowvector with simulation output timestamps.
% @param[out] Y_data: Dim x N matrix with simulation output positions.
% @param[out] dY_data: Dim x N matrix with simulation output velocities.
% @param[out] ddY_data: Dim x N matrix with simulation output accelerations.
%


if (nargin < 6)
    disturbance_fun = @(t) 0;
end

%% set initial values
if (~iscell(gmp)), gmp = {gmp}; end
Dim = length(gmp);
ddy = zeros(Dim,1);
dy = zeros(Dim,1);
y = g;
t = 0.0;
dz = zeros(Dim,1);
z = zeros(Dim,1);

t_end = abs(T);
tau = T;
x = 1.0;
x_dot = -1/tau;
s = [x; x_dot];

f_disturb = disturbance_fun(t);
F_dist = [];

iters = 0;
Time = [];
Y_data = [];
dY_data = [];
ddY_data = [];
x_data = [];

for i=1:Dim
    gmp{i}.setY0(y0(i));
    gmp{i}.setGoal(g(i));
end

%% simulate
while (true)

    %% data logging
    Time = [Time t];
    Y_data = [Y_data y];
    dY_data = [dY_data dy];  
    ddY_data = [ddY_data ddy];
    F_dist = [F_dist f_disturb];
    x_data = [x_data x];
    
    yd = zeros(Dim,1);

    %% DMP simulation
    for i=1:Dim
        y_c = 0.0;
        z_c = f_disturb;
        gmp{i}.update(s, y(i), z(i), y_c, z_c);
        dy(i) = gmp{i}.getYdot();
        dz(i) = gmp{i}.getZdot();
        
        yc_dot = 0.0;
        ddy(i) = gmp{i}.getYddot(yc_dot);
        
        yd(i) = gmp{i}.getYd(x);
        
    end
    
    f_disturb = disturbance_fun(t);
    % a_stop = -0.1*norm(y-yd);
    a_stop = -0.5*norm(f_disturb);

    %% Update phase variable
    dx = -(1/tau) * exp(a_stop);

    %% Stopping criteria
    if (x<5e-3 && norm(y-y0)<5e-2) % && norm(y-g)<5e-3 && norm(dy)<5e-3)
        break;
    end

    %% Numerical integration
    iters = iters + 1;
    t = t + dt;
    x = x + dx*dt;
    y = y + dy*dt;
    z = z + dz*dt;

    if (x < 0), x = 0; end
    s = [x; dx];
    
end


end

