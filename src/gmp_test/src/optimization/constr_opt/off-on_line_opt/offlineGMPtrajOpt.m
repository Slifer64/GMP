function [Time, P_data, dP_data, ddP_data] = offlineGMPtrajOpt(gmp0, tau, y0, yg, pos_lim, vel_lim, accel_lim, opt_pos, opt_vel, use_matlab_solver)
    
    gmp = gmp0.deepCopy();
    
    %% --------  Init sim  --------
    gmp.setScaleMethod(TrajScale_Prop(3));
    gmp.setY0(y0);
    gmp.setGoal(yg);
    
    t = 0;
    dt = 0.01;
    x = t/tau;
    x_dot = 1/tau;
    x_ddot = 0;
    
    %% --------  Init MPC  --------
    
    % control horizon
%     N = ceil(tau/dt);
    
%     n_points = 700;
%     
%     dt = tau/n_points;
    
    N = ceil(tau/dt);
    
%     dt
    N

    % stiffness and damping of the system
    % y_dddot = -K*y - D*y_dot + u
    K = 300;
    D = 80;
    
    n = 6; % state dim : x = [y y_dot]
    m = 3; % control input dim : u
    
    In = eye(n,n);
    
    % state for control horizon
    % Z = [x(1), x(2), ... x(N), u(0), u(1), ... u(N-1)]
    
    % State tracking gains: (x(i) - xd(i))'*Qi*(x(i) - xd(i))
    Qi = blkdiag( opt_pos*eye(3,3) , opt_vel*10*eye(3,3) );
    QN = 1*eye(n,n);
    % Control input minimization gain: u(i)'*Ri*u(i)
    Ri = 0*eye(m,m);
    
    H = blkdiag( kron(speye(N-1), Qi), QN, kron(speye(N), Ri) );
    
    % Discrete state transition matrices:
    % x(i+1) = Ai*x(i) + Bi*u(i)
    % x(1)             - B0*u(0) = A0*x(0)
    % x(2)   - A1*x(1) - B1*u(1) = 0
    % ...
    % x(i+1) - Ai*x(i) - Bi*u(i) = 0
    Ai = (In + [zeros(3,3) eye(3,3); -K*eye(3,3) -D*eye(3,3)]*dt);
    Bi = [zeros(3,3); eye(3,3)]*dt;
    
    % State transition equality constraints for entire horizon: 
    % Ax*Z(1:n*N) + Bu*Z(n*N+1:end) = [A0*z0, 0, 0, ..., 0]
    Ax = kron(speye(N), speye(n)) + kron(sparse(diag(ones(N-1, 1), -1)), -Ai);
    Bu = kron(speye(N), -Bi);
    % equality constraints matrix
    Aeq = [Ax, Bu];
    
    % acceleration constraints using numerical diff: x_ddot(i) = (x_dot(i+1) - x_dot(i))/dt
    D_a = sparse( [zeros(3,3), eye(3,3)/dt] );
    temp = diag(ones(N-1, 1), 1);
    Aineq_accel =  [kron(-speye(N-1,N) + sparse(temp(1:end-1,:)), D_a) , zeros((N-1)*3,m*N)];
    % accel bounds
    accel_lb = repmat(accel_lim(:,1), N-1, 1);
    accel_ub = repmat(accel_lim(:,2), N-1, 1);

    % control input bounds
    u_min = -inf;
    u_max = inf;
    z_min = [pos_lim(:,1); vel_lim(:,1)];
    z_max = [pos_lim(:,2); vel_lim(:,2)];
    
    % extend bounds for the entire horizon N
    Z_lb = [repmat(z_min, N,1); repmat(u_min*ones(m,1),N,1)];
    Z_ub = [repmat(z_max, N,1); repmat(u_max*ones(m,1),N,1)];
    
    % Add final state constraint as: x_final <= x(N) <= x_final
    Z_lb((N-1)*n+1 : n*N) = [yg; zeros(3,1)];
    Z_ub((N-1)*n+1 : n*N) = [yg; zeros(3,1)];
    
    % use this if you want to have: (u(i) - ud(i))'*Ri*(u(i) - ud(i))
    use_ud = 0;
    
    % if zero, then the OSQP solver is used
    % use_matlab_solver = 1;

    %% --------  calc Beq  --------
    % DMP phase variable
    si = x;
    si_dot = x_dot;
    si_ddot = x_ddot;

    yd_i = gmp.getYd(si);
    dyd_i = gmp.getYdDot(si, si_dot);
    ddyd_i = gmp.getYdDDot(si, si_dot, si_ddot);
    
    % Aeq*Z
    beq = zeros( n*N , 1);
    
    % J = sum_i=1:N { x(i)'*Qi*x(i) -2*(Qi*xd(i))'*x(i) }  + sum_i=0:N-1 { u(i)'*Ri*u(i) -2*(Ri*ud(i))'*u(i) }
    % qx(i*n+1 : (i+1)*n) = -Q(i+1)*xd(i+1)   (!!! because i starts from 1)
    qx = zeros(N*n,1);
    % qu(i*n+1 : (i+1)*n) = -Ri*ud(i)
    qu = zeros(N*m,1);
    
    % used as an initial guess for the solver
    % It's actually the states x(i) and control inputs u(i) that ensure
    % x(i)=xd(i) in the absesnce of the bounds constraints
    Z_prev = zeros(N*(n+m));

    for i=0:N-1

        % ideal control input so that x(i) = xd(i) in the absence of bound
        % constraints
        ud_i = K*yd_i + D*dyd_i + ddyd_i;

        beq((i*n)+1 : (i+1)*n) = use_ud * Bi*ud_i;

        qu((i*m)+1 : (i+1)*m) = use_ud *(-Ri*ud_i);

        % i := i+1, move one step forward, by integrating s(i) to get s(i+1) 
        % and in turn, get x(i+1).
        % This is because x(i) goes from 1:N, and the loop count goes 0:N-1
        % Although below we use subscript 'i' for the variables s,y,x its 
        % actually 'i+1' for the current loop 
        si = si + si_dot*dt;
        si_dot = si_dot + si_ddot*dt;
        % si_ddot = ... (if it changes too)

        if (si >= 1)
            si = 1;
            si_dot = 0;
            si_ddot = 0;
        end

        yd_i = gmp.getYd(si);
        dyd_i = gmp.getYdDot(si, si_dot);
        ddyd_i = gmp.getYdDDot(si, si_dot, si_ddot);
        
        xd_i_plus_1 = [yd_i; dyd_i];
        
        qx((i*n)+1 : (i+1)*n) = -Qi*xd_i_plus_1;
        
        Z_prev((i*n)+1 : (i+1)*n) = xd_i_plus_1;
        Z_prev(n*N+(i*m)+1 : n*N+(i+1)*m) = ud_i;

    end
    % repeat the final assignment using the final value QN
    qx((N-1)*n+1 : N*n) = -QN*[yd_i; dyd_i];
    
    % Substitute the initial value constraint 
    x0 = [y0; zeros(3,1)];
    beq(1:n) = beq(1:n) + Ai*x0;
    
    % Accumulate
    q = [qx; qu];

    %% ===========  solve optimization problem  ===========
    
    tic
    
    %% --------- OSQP solver ----------
    if (~use_matlab_solver)

        % Create an OSQP object
        prob = osqp;

        % - OSQP constraints
        % lb <= I*Z <= ub
        Aineq = speye(N*(n+m));
        
        % Transform equality constraints in inequalities:
        % Aeq*Z = beq : beq <= Aeq*Z <= beq
        A = [Aeq; Aineq];
        lb = [beq; Z_lb];
        ub = [beq; Z_ub];

        % Setup workspace
        prob.setup(H, q, A, lb, ub, 'warm_start',true, 'verbose',false, 'eps_abs',1e-6, 'eps_rel',1e-6);

        res = prob.solve();

        if ( res.info.status_val ~= 1)
            res.info
            if (abs(res.info.status_val) == 2), warning(res.info.status);
            else, error(res.info.status);
            end
        end
        
        X = reshape( res.x(1:N*n), n, N );
        U = reshape( res.x(N*n+1:end), m, N );
        
    end

    %% --------- matlab solver ----------
    if (use_matlab_solver)

        A = [Aineq_accel; -Aineq_accel];
        b = [accel_ub; -accel_lb];
        
        solver_opt = optimoptions('quadprog', 'Algorithm','interior-point-convex', 'StepTolerance',0, 'Display','off', 'MaxIterations',2000);
        [Z, ~, ex_flag, opt_output] = quadprog(H,q, A,b, Aeq,beq, Z_lb,Z_ub, Z_prev, solver_opt);

        if (ex_flag == 1 || ex_flag == 2)
            % success
        else
            warning(opt_output.message);
        end
        
        X = reshape( Z(1:N*n), n, N );
        U = reshape( Z(N*n+1:end), m, N);
    end
    
    fprintf('===> Optimization finished! Elaps time: %f ms\n',toc()*1000);
        
    X = [x0 X];
   
    %% --------  Log data  --------
    Time = (0:N)*dt;
    P_data = X(1:3,:);
    dP_data = X(4:6,:);
    ddP_data = zeros(size(dP_data));
    for i=1:3
        ddP_data(i,:) = diff( [dP_data(i,:) dP_data(i,end)] ) / dt;
    end
    
%     % This should produce the same output with the above
%     X = zeros(n,N+1);
%     X(:,1) = x0;
%     xi = x0;
%     for i=1:N
%         xi = Ai*xi + Bi*U(:,i);
%         X(:,i+1) = xi;
%     end
%     P_data = X(1:3,:);
%     dP_data = X(4:6,:);
%     ddP_data = zeros(size(dP_data));
%     for i=1:3
%         ddP_data(i,:) = diff( [dP_data(i,:) dP_data(i,end)] ) / dt;
%     end
    

%     u_norm = arrayfun(@(j) norm(U(:,j)) , 1:size(U,2));
%     figure;
%     plot(Time(1:end-1), u_norm, 'LineWidth',2, 'Color','magenta');
    
%     gmp2 = GMP(3, 100, 1.5);
%     gmp2.train('LS', Time/Time(end), P_data);
%     [Time, P_data, dP_data, ddP_data] = getGMPTrajectory(gmp2, tau, y0, yg);
    

end
