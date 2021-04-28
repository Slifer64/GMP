function [Time, Y_data, dY_data, ddY_data] = simulateGMP(gmp, y0, g, T, dt)
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


%% set initial values
if (~iscell(gmp)), gmp = {gmp}; end
Dim = length(gmp);
y = y0;
dy = zeros(Dim,1);
ddy = zeros(Dim,1);
z = zeros(Dim,1);
dz = zeros(Dim,1);
t = 0.0;

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
    % x_data = [x_data x];

    %% DMP simulation
    for i=1:Dim
        y_c = 0.0;
        z_c = 0.0;
        gmp{i}.update(s, y(i), z(i), y_c, z_c);
        dy(i) = gmp{i}.getYdot();
        dz(i) = gmp{i}.getZdot();
        
        yc_dot = 0.0;
        ddy(i) = gmp{i}.getYddot(yc_dot);
    end

    %% Update phase variable
    x_dot = 1/tau;

    %% Stopping criteria
    if (t>=1.1*t_end) % && norm(y-g)<5e-3 && norm(dy)<5e-3)
        break;
    end

    %% Numerical integration
    iters = iters + 1;
    t = t + dt;
    x = x + x_dot*dt;
    y = y + dy*dt;
    z = z + dz*dt;

    s = [x; x_dot];
    
end


end

