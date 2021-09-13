function [Time, P_data, dP_data, ddP_data] = onlineGMPTrajOpt(gmp, tau, y0, yg, pos_lim, vel_lim, accel_lim)
    
    gmp2 = gmp.deepCopy();
    
    %% --------  Init sim  --------
    gmp2.setScaleMethod(TrajScale_Prop(3));
    gmp2.setY0(y0);
    gmp2.setGoal(yg);

    Time = [];
    P_data = [];
    dP_data = [];
    ddP_data = [];

    p = y0;
    p_dot = zeros(size(p));
    p_ddot = zeros(size(p));

    t = 0;
    dt = 0.002;
    x = t/tau;
    x_dot = 1/tau;
    x_ddot = 0;
    
    %% --------  Init MPC  --------
    N = 30;
    
    K = 300;
    D = 80;
    
    n = 6; % state dim
    m = 3; % control input dim
    
    In = eye(n,n);
    
    % H = eye((n+m)*N, (n+m)*N);
    Qi = blkdiag( 1*eye(3,3) , 1*eye(3,3) );
    QN = 100*eye(n,n);
    Ri = 0*eye(m,m);

    %H = kron(eye(N,N), blkdiag(Qi,Ri) ); % speye?
    H = blkdiag( kron(speye(N-1), Qi), QN, kron(speye(N), Ri) );
    
    Ai = (In + [zeros(3,3) eye(3,3); -K*eye(3,3) -D*eye(3,3)]*dt);
    Bi = [zeros(3,3); eye(3,3)]*dt;
    
    Ax = kron(speye(N), speye(n)) + kron(sparse(diag(ones(N-1, 1), -1)), -Ai);
    Bu = kron(speye(N), -Bi);
    Aeq = [Ax, Bu];

    u_min = -inf;
    u_max = inf;
    z_min = [pos_lim(:,1); vel_lim(:,1)];
    z_max = [pos_lim(:,2); vel_lim(:,2)];
    
    X_lb = [repmat(z_min, N,1); repmat(u_min*ones(m,1),N,1)];
    X_ub = [repmat(z_max, N,1); repmat(u_max*ones(m,1),N,1)];
    
    Aineq = speye(N*(n+m));
    
    solver_opt = optimoptions('quadprog', 'Algorithm','interior-point-convex', 'StepTolerance',0, 'Display','off', 'MaxIterations',2000);
    X_prev = zeros(N*(n+m),1);
    
    % Create an OSQP object
    prob = osqp;
    
    % - OSQP constraints
    q = zeros(N*(n+m),1);
    beq = zeros( n*N , 1);
    
    A = [Aeq; Aineq];
    lb = [beq; X_lb];
    ub = [beq; X_ub];

    % Setup workspace
    prob.setup(H, q, A, lb, ub, 'warm_start',true, 'verbose',false, 'eps_abs',1e-6, 'eps_rel',1e-6);
    
    use_ud = 0;
    
    use_matlab_solver = 0;

    %% ===========  Simulation loop  ===========
    while (true)
        
        %% --------  Stopping criteria  --------
        if (x > 1.0), break; end
        
        t/tau
        
        if (x >= 1)
            x_dot = 0;
        end

        
        %% --------  calc Beq  --------
        qx = zeros(N*n,1);
        qu = zeros(N*m,1);

        xi = x;
        xi_dot = x_dot;
        xi_ddot = x_ddot;
        
        beq = zeros( n*N , 1);
            
        pd_i = gmp2.getYd(xi);
        dpd_i = gmp2.getYdDot(xi, xi_dot);
        ddpd_i = gmp2.getYdDDot(xi, xi_dot, xi_ddot);
            
        j = n*N + 1;
        
        xi_one_reached = false;
        
        for i=0:N-1

            ud_i = K*pd_i + D*dpd_i + ddpd_i;
            
            beq((i*n)+1 : (i+1)*n) = use_ud * Bi*ud_i;
            
            qu((i*m)+1 : (i+1)*m) = zeros(m,1); % -Ri*ud_i;
            
            % i := i+1
            xi = xi + xi_dot*dt;
            xi_dot = xi_dot + xi_ddot*dt;
            
            if (xi >= 1)
                xi = 1;
                xi_dot = 0;
                xi_ddot = 0;
            end
            
            pd_i = gmp2.getYd(xi);
            dpd_i = gmp2.getYdDot(xi, xi_dot);
            ddpd_i = gmp2.getYdDDot(xi, xi_dot, xi_ddot);
            
            if (xi >= 1) % && ~xi_one_reached) % we reached the target
                %norm([pd_i; dpd_i] - [yg; zeros(3,1)])
%                 lb(j:j+n-1,:) = [pd_i; dpd_i] - 2e-2; %[yg; zeros(3,1)];
%                 ub(j:j+n-1,:) = [pd_i; dpd_i] + 2e-2;
%                 
%                 X_lb(i*n+1 : (i+1)*n) = [pd_i; dpd_i] - 16e-2;
%                 X_ub(i*n+1 : (i+1)*n) = [pd_i; dpd_i] + 16e-2;
                
                xi_one_reached = true;
            end
            j = j + n;
            
            qx((i*n)+1 : (i+1)*n) = -Qi*[pd_i; dpd_i];

        end
        qx((N-1)*n+1 : N*n) = -QN*[pd_i; dpd_i];
        z0 = [p; p_dot];
        beq(1:n) = beq(1:n) + Ai*z0;
        
        q = [qx; qu];

        lb(1:n*N) = beq;
        ub(1:n*N) = beq;
        
        %% --------  calc auxiliary control u  --------
        
        u = 0;
        
        if (~use_matlab_solver)

            prob.update('q', q, 'l', lb, 'u', ub);

            res = prob.solve();

            if ( res.info.status_val ~= 1)
                res.info
                if (abs(res.info.status_val) == 2), warning(res.info.status);
                else, error(res.info.status);
                end
            end

            u = res.x(N*n+1:N*n+m);
        end
        
%         if (~ isempty( find( Aineq * res.x > X_ub +1e-3 ) ) )
%             ind = find( Aineq * res.x < X_ub );
%             temp = Aineq * res.x - X_lb;
%             temp(temp>0)
%             res.info
%             warning('X_ub Constraints violation!'); 
%         end
%         
%         if (~ isempty( find( Aineq * res.x < X_lb -1e-3 ) ) )
%             temp = Aineq * res.x - X_lb;
%             temp(temp<0)
%             res.info
%             warning('X_lb Constraints violation!'); 
%         end
        
        if (use_matlab_solver)
            [X, ~, ex_flag, opt_output] = quadprog(H,q, [],[], Aeq,beq, X_lb,X_ub, X_prev, solver_opt);

            if (ex_flag == 1 || ex_flag == 2)
                % success
            else
                warning(opt_output.message);
            end

            u = X(N*n+1:N*n+m);

            X_prev = [X(n+m+1:end); X(end-m-n+1:end)];
        end

        %% --------  calc u_ref  --------
        p_ref = gmp2.getYd(x);
        dp_ref = gmp2.getYdDot(x, x_dot);
        ddp_ref = gmp2.getYdDDot(x, x_dot, x_ddot);
        
        u_ref = use_ud * ( K*p_ref + D*dp_ref + ddp_ref);
        
        %% --------  Simulate dynamics  --------
        p_ddot = -K*p - D*p_dot + u_ref + u;
        
        z = [p; p_dot];
        z2 = Ai*z + Bi*(u_ref + u);
        
        p_ddot = ( z2(4:6) - z(4:6) ) / dt;
        
        %% --------  Log data  --------
        Time = [Time t];
        P_data = [P_data p];
        dP_data = [dP_data p_dot];
        ddP_data = [ddP_data p_ddot];
    
        %% --------  Numerical integration  --------
        t = t + dt;
        x = x + x_dot*dt;
        x_dot = x_dot + x_ddot*dt;
%         p = p + p_dot*dt;
%         p_dot = p_dot + p_ddot*dt;
        p = z2(1:3);
        p_dot = z2(4:6);

    end
    
    p_err = abs(p - yg)'
    dp_err = abs(p_dot - zeros(3,1))'
        
end
