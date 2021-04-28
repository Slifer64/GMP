function [Time, Q_data, rotVel_data, rotAccel_data] = simulateGMPo_in_quat_space(gmp_o, Q0, Qg, T, dt)
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

%% simulate
while (true)

    %% data logging
    Time = [Time t];
    Q_data = [Q_data Q];
    rotVel_data = [rotVel_data rotVel];  
    rotAccel_data = [rotAccel_data rotAccel];

    
    %% GMP simulation
    yc = 0; % optional coupling for 'y' state
    zc = 0; % optional coupling for 'z' state
    yc_dot = 0; % derivative of coupling for 'y' state
    s = [x; x_dot; x_ddot];
    rotAccel = gmp_o.calcRotAccel(s, Q, rotVel, Qg, yc, zc, yc_dot);

    
    %% Stopping criteria   
    if (t>1.5*t_end)
        warning('Time limit reached... Stopping simulation!');
        break;
    end
    
    eo = quatLog(quatProd(Qg, quatInv(Q)));
    if (t>=t_end && norm(eo)<0.02)
        break;
    end

    
    %% Numerical integration
    t = t + dt;
    x = x + x_dot*dt;
    Q = quatProd( quatExp(rotVel*dt), Q);
    rotVel = rotVel + rotAccel*dt;  
    
end

end

