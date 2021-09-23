function [Time, P_data, dP_data, ddP_data] = onlineGMPtrajOpt(gmp0, tau, y0, yg, pos_lim, vel_lim, accel_lim, opt_pos, opt_vel, qp_solver_type)
    
    gmp = gmp0.deepCopy();
    
    n_dof = length(y0);
    
    %% --------  Init sim  --------
    gmp.setScaleMethod(TrajScale_Prop(n_dof));
    gmp.setY0(y0);
    gmp.setGoal(yg);
    
    Time = [];
    P_data = [];
    dP_data = [];
    ddP_data = [];

    y = y0;
    y_dot = zeros(size(y));
    y_ddot = zeros(size(y));

    t = 0;
    dt = 0.005;
    s = t/tau;
    s_dot = 1/tau;
    s_ddot = 0;
    
    x_final = [yg; zeros(n_dof,1)];
    
    %% --------  Init MPC  --------
    N = 10;%10; %200;
    dt_ = dt;%dt;%0.015;%dt; %0.010;
    use_ud = 0;

    % stiffness and damping of the system
    % y_dddot = -K*y - D*y_dot + u
    K = 300;
    D = 80;
    
    n = 2*n_dof; % state dim : x = [y y_dot]
    m = n_dof; % control input dim : u
    
    In = eye(n,n);
    
    % state for control horizon
    % Z = [x(1), x(2), ... x(N), u(0), u(1), ... u(N-1)]
    
    % State tracking gains: (x(i) - xd(i))'*Qi*(x(i) - xd(i))
    Qi = blkdiag( opt_pos*eye(n_dof,n_dof) , opt_vel*10*eye(n_dof,n_dof) );
    QN = Qi; %1*eye(n,n);
    % Control input minimization gain: u(i)'*Ri*u(i)
    Ri = 0*eye(m,m);
    
    H = blkdiag( kron(speye(N-1), Qi), QN, kron(speye(N), Ri) );
    %H = 1e-12*speye(size(H)) + H;
    
    % Discrete state transition matrices:
    % x(i+1) = Ai*x(i) + Bi*u(i)
    % x(1)             - B0*u(0) = A0*x(0)
    % x(2)   - A1*x(1) - B1*u(1) = 0
    % ...
    % x(i+1) - Ai*x(i) - Bi*u(i) = 0
    
    % keep the simulation time-step 'dt' for the first prediction
    A0 = (In + [zeros(n_dof,n_dof) eye(n_dof,n_dof); -K*eye(n_dof,n_dof) -D*eye(n_dof,n_dof)]*dt);
    B0 = [zeros(n_dof,n_dof); eye(n_dof,n_dof)]*dt;
    % use a different time-step 'dt_' for the next predictions
    Ai = (In + [zeros(n_dof,n_dof) eye(n_dof,n_dof); -K*eye(n_dof,n_dof) -D*eye(n_dof,n_dof)]*dt_);
    Bi = [zeros(n_dof,n_dof); eye(n_dof,n_dof)]*dt_;
    
    % State transition equality constraints for entire horizon: 
    % Ax*Z(1:n*N) + Bu*Z(n*N+1:end) = [A0*x0, 0, 0, ..., 0]
    Ax = kron(speye(N), speye(n)) + kron(sparse(diag(ones(N-1, 1), -1)), -Ai);
    Bu = blkdiag(-B0, kron(speye(N-1), -Bi) );
    % equality constraints matrix
    Aeq = [Ax, Bu];
    
    % acceleration constraints using numerical diff: x_ddot(i) = (x_dot(i+1) - x_dot(i))/dt
    D_a = sparse( [zeros(n_dof,n_dof), eye(n_dof,n_dof)] );
    temp = diag(ones(N-1, 1), 1);
    Aineq_accel =  [zeros(n_dof), eye(n_dof) zeros(n_dof, N*(m+n)-n);
        kron(-speye(N-1,N) + sparse(temp(1:end-1,:)), D_a) , zeros((N-1)*n_dof,m*N)];
    % accel bounds
    % constraint for initial acceleration: a_lb*dt <= y_dot(1) - y_dot(0) <= a_ub*dt
    y0_dot = y_dot;
    accel_lb = [accel_lim(:,1)*dt + y0_dot; repmat(accel_lim(:,1)*dt_, N-1, 1)];
    accel_ub = [accel_lim(:,2)*dt + y0_dot; repmat(accel_lim(:,2)*dt_, N-1, 1)];

    % control input bounds
    u_min = -inf(m,1);
    u_max = inf(m,1);
    x_min = [pos_lim(:,1); vel_lim(:,1)];
    x_max = [pos_lim(:,2); vel_lim(:,2)];
    
    % extend bounds for the entire horizon N
    Z_lb = [repmat(x_min, N,1); repmat(u_min,N,1)];
    Z_ub = [repmat(x_max, N,1); repmat(u_max,N,1)];
    
%     % Add final state constraint as: x_final <= x(N) <= x_final
%     Z_lb((N-1)*n+1 : n*N) = [yg; zeros(n_dof,1)];
%     Z_ub((N-1)*n+1 : n*N) = [yg; zeros(n_dof,1)];

    %% --------  Init solver  --------
    if (qp_solver_type == 1)
        
        % Create an OSQP object
        prob = osqp;

        % - OSQP constraints
        % lb <= I*Z <= ub
        Aineq = speye(N*(n+m));

        % Transform equality constraints in inequalities:
        % Aeq*Z = beq : beq <= Aeq*Z <= beq
        A_osqp = [Aeq; Aineq; Aineq_accel];
        
        osqp_initialized = false;
        
    else
        
        solver_opt = optimoptions('quadprog', 'Algorithm','interior-point-convex', 'StepTolerance',0, 'Display','off', 'MaxIterations',2000);
        
        % used as an initial guess for the solver
        % It's actually the states x(i) and control inputs u(i) that ensure
        % x(i)=xd(i) in the absesnce of the bounds constraints
        Z_init = zeros(N*(n+m));

        si = s;
        si_dot = s_dot;
        si_ddot = s_ddot;

        yd_i = gmp.getYd(si);
        dyd_i = gmp.getYdDot(si, si_dot);
        ddyd_i = gmp.getYdDDot(si, si_dot, si_ddot);
            
        si = si + si_dot*dt;
        si_dot = si_dot + si_ddot*dt;
            
        for i=1:N

            % ideal control input so that x(i) = xd(i) in the absence of bound
            % constraints
            ud_i_minus_1 = K*yd_i + D*dyd_i + ddyd_i;

            if (si >= 1)
                si = 1;
                si_dot = 0;
                si_ddot = 0;
            end

            yd_i = gmp.getYd(si);
            dyd_i = gmp.getYdDot(si, si_dot);
            ddyd_i = gmp.getYdDDot(si, si_dot, si_ddot);

            xd_i = [yd_i; dyd_i];

            Z_init((i-1)*n+1 : i*n) = xd_i;
            Z_init(n*N+(i-1)*m+1 : n*N+i*m) = ud_i_minus_1;
            
            si = si + si_dot*dt_;
            si_dot = si_dot + si_ddot*dt_;

        end
        
        
    end
    
    text_prog = ProgressText(40);
    text_prog.init();
    
    t_start = tic;
    
    u_data = [];

    %% --------  Simulation loop  --------
    while (true)
        
        %% --------  Stopping criteria  --------
        if (s > 1.0), break; end
        
        %t/tau
        text_prog.update(100*t/tau);
%         fprintf('progress: %.1f\n',t/tau);
        
        if (s >= 1)
            s_dot = 0;
        end

        %% -------  calc Beq  --------
        % DMP phase variable
        si = s;
        si_dot = s_dot;
        si_ddot = s_ddot;

        yd_i = gmp.getYd(si);
        dyd_i = gmp.getYdDot(si, si_dot);
        ddyd_i = gmp.getYdDDot(si, si_dot, si_ddot);
        
        ud = K*yd_i + D*dyd_i + ddyd_i;

        % Aeq*Z
        beq = zeros( n*N , 1);

        % J = sum_i=1:N { x(i)'*Qi*x(i) -2*(Qi*xd(i))'*x(i) }  + sum_i=0:N-1 { u(i)'*Ri*u(i) -2*(Ri*ud(i))'*u(i) }
        % qx(i*n+1 : (i+1)*n) = -Q(i+1)*xd(i+1)   (!!! because i starts from 1)
        qx = zeros(N*n,1);
        % qu(i*n+1 : (i+1)*n) = -Ri*ud(i)
        qu = zeros(N*m,1);
        
        % restore bounds in case they were changed in previous iter for si>=1
        Z_lb(1:n*N) = repmat(x_min, N,1);
        Z_ub(1:n*N) = repmat(x_max, N,1);
        
        fc_flag = 0;
        
        % move one step forward to get the next desired state
        si = si + si_dot*dt;
        si_dot = si_dot + si_ddot*dt;

        for i=1:N

            if (si >= 1)
                si = 1;
                si_dot = 0;
                si_ddot = 0;
                %if (i==N-1)
                if (~fc_flag)
                    Z_lb((i-1)*n+1:i*n) = x_final;
                    Z_ub((i-1)*n+1:i*n) = x_final;
                    fc_flag = 1;
                end
                
                yd_i = yg;
                dyd_i = zeros(n_dof,1);
                ddyd_i = zeros(n_dof,1);

            else
                yd_i = gmp.getYd(si);
                dyd_i = gmp.getYdDot(si, si_dot);
                ddyd_i = gmp.getYdDDot(si, si_dot, si_ddot);
            end

            xd_i = [yd_i; dyd_i];
            qx( (i-1)*n + 1 : i*n ) = -Qi*xd_i;
            
            % qu((i-1)*m+1 : i*m) = -Ri*ud_i_minus_1;
            
            if (i<N)
                ud_i = K*yd_i + D*dyd_i + ddyd_i;
                beq(i*n+1 : (i+1)*n) = use_ud*Bi*ud_i;
            end
            
            % move one step forward
            si = si + si_dot*dt_;
            si_dot = si_dot + si_ddot*dt_;
        end
        % repeat the final assignment using the final value QN
        qx((N-1)*n+1 : N*n) = -QN*[yd_i; dyd_i];

        % Substitute the initial value constraint 
        x0 = [y; y_dot];
        beq(1:n) = A0*x0 + use_ud*B0*ud;
        
        % update initial acceleration constraint
        y0_dot = x0(n_dof+1:end);
        accel_lb(1:n_dof) = accel_lim(:,1)*dt + y0_dot;
        accel_ub(1:n_dof) = accel_lim(:,2)*dt + y0_dot;

        % accumulate
        q = [qx; qu];

        %% ===========  solve optimization problem  ==========

        %% --------- OSQP solver ----------
        if (qp_solver_type == 1)
            
            lb = [beq; Z_lb; accel_lb];
            ub = [beq; Z_ub; accel_ub];

            if (~osqp_initialized)
                % Setup workspace
                prob.setup(H, q, A_osqp, lb, ub, 'warm_start',true, 'verbose',false); %, 'eps_abs',1e-6, 'eps_rel',1e-6, 'max_iter',30000);
                osqp_initialized = true;
            else
                prob.update('q', q, 'l', lb, 'u', ub);
            end
            
            res = prob.solve();

            if ( res.info.status_val ~= 1)
                res.info
                if (abs(res.info.status_val) == 2), warning(res.info.status);
                else, error(res.info.status);
                end
                text_prog.printInNewLine();
            end
            
            u = res.x(N*n+1:N*n+m);

        end

        %% --------- matlab solver ----------
        if (qp_solver_type == 0)

            A = [Aineq_accel];%; -Aineq_accel];
            b = [accel_ub];%; -accel_lb];

            [Z, ~, ex_flag, opt_output] = quadprog(H,q, A,b, Aeq,beq, Z_lb,Z_ub, Z_init, solver_opt);

            if (ex_flag == 1 || ex_flag == 2)
                % success
            else
                warning(opt_output.message);
                text_prog.printInNewLine();
                if (ex_flag < 0), break; end
            end
            
            u = Z(N*n+1:N*n+m);
            
            % use the previous solution as initial guess for next iter
            %Z_init = Z;
            
            % discard x(1), u(0) and replicate x(N) and u(N-1)
            Z_init = [Z(n+1:n*N); Z((N-1)*n+1:n*N); Z(n*N+m+1:end); Z(end-m+1:end)]; 

        end
        
        %% --------  Simulate dynamics  --------
        %y_ddot = -K*y - D*y_dot + u;
        
        x = [y; y_dot];
        x_next = A0*x + B0*(u + use_ud*ud);
        
        u_data = [u_data norm(u)];
        
        y_ddot = ( x_next(n_dof+1:n) - x(n_dof+1:n) ) / dt;
        
        %% --------  Log data  --------
        Time = [Time t];
        P_data = [P_data y];
        dP_data = [dP_data y_dot];
        ddP_data = [ddP_data y_ddot];
    
        %% --------  Numerical integration  --------
        t = t + dt;
        s = s + s_dot*dt;
        s_dot = s_dot + s_ddot*dt;
%         y = y + y_dot*dt;
%         y_dot = y_dot + y_ddot*dt;
        y = x_next(1:n_dof);
        y_dot = x_next(n_dof+1:n);
    
    end
    
%     ddP_data = zeros(n_dof,length(Time));
%     for i=1:n_dof, ddP_data(i,:) = [diff(dP_data(i,:)) 0]/dt; end
    
    text_prog.update(100);
    fprintf('\n');
    fprintf('===> online-GMP-traj optimization finished! Elaps time: %f ms\n',toc(t_start)*1000);
    
    
    figure;
    plot(Time, u_data, 'LineWidth',2, 'Color','red')

end
