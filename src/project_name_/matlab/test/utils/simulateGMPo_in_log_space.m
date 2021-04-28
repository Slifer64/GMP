function [Time, Q_data, rotVel_data, rotAccel_data] = simulateGMPo_in_log_space(gmp_o, Q0, Qg, T, dt)
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
Q_prev = Q;
rotVel = zeros(3,1);
rotAccel = zeros(3,1);
q = gmp_o.quat2q(Q0, Q0);
qdot = zeros(3,1);
dy = zeros(3,1);
dz = zeros(3,1);

gmp_o.setQ0(Q0);
gmp_o.setQg(Qg);
y = gmp_o.getY(Q);
z = gmp_o.getZ(rotVel, Q);
g = gmp_o.quat2q(Qg, Q0);


%% simulate
while (true)

    %% data logging
    Time = [Time t];
    Q_data = [Q_data Q];
    rotVel_data = [rotVel_data rotVel];  
    rotAccel_data = [rotAccel_data rotAccel];
    
    yc_dot = 0;

    %% DMP simulation
    yc = zeros(3,1);
    zc = zeros(3,1);
    s = [x; x_dot; x_ddot];
    gmp_o.update(s, y, z, yc, zc);

    dy = gmp_o.getYdot();
    dz = gmp_o.getZdot();
    rotAccel = gmp_o.getRotAccel(Q, yc_dot);
    

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
    y = y + dy*dt;
    z = z + dz*dt;
    
    q = y;
    dy = z;
    qdot = dy;
    
    Q_prev = Q;
    Q = gmp_o.q2quat(q, Q0);
    if (Q_prev'*Q<0), Q = -Q; end
    
    Q1 = gmp_o.quatTf(Q, Q0);
    rotVel = gmp_o.qdot2rotVel(qdot, Q1);
    
end


end

