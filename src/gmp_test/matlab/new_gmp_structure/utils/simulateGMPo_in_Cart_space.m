function [Time, Q_data, rotVel_data, rotAccel_data] = simulateGMPo_in_Cart_space(gmp_o, Q0, Qg, T, dt)
%% Simulates a dmp encoding Cartesian orientation usning unit quaternions.


%% set initial values
t_end = T;
tau = t_end;

Time = [];
Q_data = [];
rotVel_data = [];
rotAccel_data = [];

t = 0.0;
x = 0.0;
x_dot = 1/tau;
x_ddot = 0;
Q = Q0;
rotVel = zeros(3,1);
rotAccel = zeros(3,1);

gmp_o.setQ0(Q0);   % set initial orientation
gmp_o.setQg(Qg);   % set target orientation

elaps_t = [];

%% simulate
while (true)

    %% data logging
    Time = [Time t];
    Q_data = [Q_data Q];
    rotVel_data = [rotVel_data rotVel];  
    rotAccel_data = [rotAccel_data rotAccel];

    tic;
    
    %% GMP simulation
%     Qd = gmp_o.getQd(x);
%     Vd = gmp_o.getVd(x, x_dot);
%     Vd_dot = gmp_o.getVdDot(x, x_dot, x_ddot);
    [Qd, Vd, Vd_dot] = gmp_o.getRefTraj(x, x_dot, x_ddot);
    
    rotAccel = Vd_dot + 5*(Vd-rotVel) + 20*math_.quatLog(math_.quatDiff(Qd,Q));
    
    elaps_t = [elaps_t toc()*1000];
    
    %% Stopping criteria   
    if (t>1.5*t_end)
        warning('Time limit reached... Stopping simulation!');
        break;
    end
    
    eo = gmp_.quatLog(gmp_.quatProd(Qg, gmp_.quatInv(Q)));
    if (t>=t_end && norm(eo)<0.02)
        break;
    end

    %% Numerical integration
    t = t + dt;
    x = x + x_dot*dt;
    Q = gmp_.quatProd( gmp_.quatExp(rotVel*dt), Q);
    rotVel = rotVel + rotAccel*dt;  
    
end

mean_elaps_t = mean(elaps_t)
std_elaps_t = std(elaps_t,1)

end

