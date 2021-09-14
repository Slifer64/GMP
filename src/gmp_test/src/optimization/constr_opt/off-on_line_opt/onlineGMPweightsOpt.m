function [Time, P_data, dP_data, ddP_data] = onlineGMPweightsOpt(gmp0, tau, y0, yg, pos_lim, vel_lim, accel_lim, opt_pos, opt_vel, use_matlab_solver)
    
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
    
    x_final = [yg; zeros(n_dof,1); zeros(n_dof,1)];
    n_dof3 = length(x_final); % for pos, vel, accel
    
    %% --------  Init MPC  --------
    N = 8;%10; %200;
    
    
    N_kernels = gmp.numOfKernels();
    
    n = n_dof * N_kernels;
    
    % state for control horizon
    % Z = [x(1), x(2), ... x(N), u(0), u(1), ... u(N-1)]
    
    % State tracking gains: (x(i) - xd(i))'*Qi*(x(i) - xd(i))
    Qi = blkdiag( opt_pos*eye(n_dof,n_dof) , opt_vel*10*eye(n_dof,n_dof) );
    QN = Qi; %blkdiag( 10*eye(n_dof,n_dof) , 0*eye(n_dof,n_dof) );

    z_min = [pos_lim(:,1); vel_lim(:,1); accel_lim(:,1)];
    z_max = [pos_lim(:,2); vel_lim(:,2); accel_lim(:,2)];
    
    w = gmp.W';
    w = w(:); % initial guess for weights

    %% --------  Init solver  --------
    if (~use_matlab_solver)
        
        % Create an OSQP object
        prob = osqp;
        
        osqp_initialized = false;
        
    else
        
        solver_opt = optimoptions('quadprog', 'Algorithm','interior-point-convex', 'StepTolerance',0, 'Display','off', 'MaxIterations',2000);

    end
    
    text_prog = ProgressText(40);
    text_prog.init();
    
    tic

    y0_ = y0;
    y0_dot = zeros(n_dof,1);
    y0_ddot = zeros(n_dof,1);
    phi0 = gmp.regressVec(s);
    phi0_dot = gmp.regressVecDot(s, s_dot);
    phi0_ddot = gmp.regressVecDDot(s, s_dot, s_ddot);
    
    phi_f = gmp.regressVec(1);
    phi_f_dot = gmp.regressVecDot(1, 0);
    phi_f_ddot = gmp.regressVecDDot(1, 0, 0);

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

        %% -------  calc problem matrices  --------
        
%         % restore bounds in case they were changed in previous iter for si>=1
%         Z_min = repmat(z_min, N,1);
%         Z_max = repmat(z_max, N,1);
        
        H = sparse(n,n);
        q = zeros(n,1);
        Aineq = zeros(N*n_dof3, n);
        
        % DMP phase variable
        si = s;
        si_dot = s_dot;
        si_ddot = s_ddot;
        dt_ = 1.5*dt; %0.8*dt;
        
        for i=1:N
            
            if (si >= 1)
                N = i-1;
                break;
            end

            yd_i = gmp.getYd(si);
            dyd_i = gmp.getYdDot(si, si_dot);
            
            phi = gmp.regressVec(si);
            phi_dot = gmp.regressVecDot(si, si_dot);
            phi_ddot = gmp.regressVecDDot(si, si_dot, si_ddot);

            Psi = [kron(eye(n_dof),phi'); kron(eye(n_dof),phi_dot')];
            
            if (i==N), Qi_ = QN;
            else, Qi_ = Qi;
            end
            
            H = H + Psi'*Qi_*Psi;
            q = q - Psi'*Qi_*[yd_i; dyd_i];
            
            Aineq((i-1)*n_dof3+1 : i*n_dof3, :) = [kron(eye(n_dof),phi'); kron(eye(n_dof),phi_dot'); kron(eye(n_dof),phi_ddot')];
            
            si = si + si_dot*dt_;
            si_dot = si_dot + si_ddot*dt_;
            % si_ddot = ... (if it changes too)
        end
        
        H = (H+H')/2; % to account for numerical errors
        
        
        Z_min = repmat(z_min, N,1);
        Z_max = repmat(z_max, N,1);
        
        Aineq = sparse(Aineq(1:N*n_dof3,:));
        
        Psi0 = [kron(eye(n_dof),phi0')];%; kron(eye(n_dof),phi0_dot')];%; kron(eye(n_dof),phi0_ddot')];
        Aeq = sparse(Psi0);
        beq = [y0_];%; y0_dot];%; y0_ddot];
        
        if (si >= 1)
            
            Phi_f = [kron(eye(n_dof),phi_f')];%; kron(eye(n_dof),phi_f_dot'); kron(eye(n_dof),phi_f_ddot')];
            Aeq = [Aeq; Phi_f];
            beq = [beq; x_final(1:n_dof)];
            
        end
        
        %% ===========  solve optimization problem  ==========

        %% --------- OSQP solver ----------
        if (~use_matlab_solver)
            
            lb = [beq; Z_min];
            ub = [beq; Z_max];

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
        if (use_matlab_solver)

            A = [Aineq; -Aineq];
            b = [Z_max; -Z_min];

            [w, ~, ex_flag, opt_output] = quadprog(H,q, A,b, Aeq,beq, [],[], w, solver_opt);

            if (ex_flag == 1 || ex_flag == 2)
                % success
            else
                warning(opt_output.message);
                text_prog.printInNewLine();
                if (ex_flag < 0), break; end
            end
            
            W = reshape(w, N_kernels, n_dof)';

        end
        
        %% --------  Simulate dynamics  --------
        
        y = W*gmp.regressVec(s);
        y_dot = W*gmp.regressVecDot(s, s_dot);
        y_ddot = W*gmp.regressVecDDot(s, s_dot, s_ddot);
        
        y0_ = y;
        y0_dot = y_dot;
        y0_ddot = y_ddot;
        phi0 = gmp.regressVec(s);
        phi0_dot = gmp.regressVecDot(s, s_dot);
        phi0_ddot = gmp.regressVecDDot(s, s_dot, s_ddot);
        
        %% --------  Log data  --------
        Time = [Time t];
        P_data = [P_data y];
        dP_data = [dP_data y_dot];
        ddP_data = [ddP_data y_ddot];
    
        %% --------  Numerical integration  --------
        t = t + dt;
        s = s + s_dot*dt;
        s_dot = s_dot + s_ddot*dt;
    
    end
    
    text_prog.update(100);
    fprintf('\n');
    fprintf('===> online-GMP-weights optimization finished! Elaps time: %f ms\n',toc()*1000);

end
