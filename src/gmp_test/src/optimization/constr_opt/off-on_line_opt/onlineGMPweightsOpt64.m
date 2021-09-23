function [Time, P_data, dP_data, ddP_data] = onlineGMPweightsOpt(gmp0, tau, y0, yg, pos_lim, vel_lim, accel_lim, opt_pos, opt_vel, use_matlab_solver)
    
    gmp = gmp0.deepCopy();
    
    n_dof = length(y0);
    
    %% --------  Init sim  --------
    gmp.setScaleMethod(TrajScale_Prop(n_dof));
    gmp.setY0(y0);
    gmp.setGoal(yg);
    
    N_kernels = 10; %gmp0.numOfKernels();
    s_data = 0:0.01:1;
    Yd_data = zeros(n_dof, length(s_data));
    for j=1:length(s_data), Yd_data(:,j)=gmp.getYd(s_data(j)); end
    gmp = GMP(n_dof, N_kernels, 1.5);
    gmp.train('LS', s_data, Yd_data);
    
    Time = [];
    P_data = [];
    dP_data = [];
    ddP_data = [];

    y = y0;
    y_dot = zeros(size(y));
    y_ddot = zeros(size(y));
    
    K = 300;
    D = 80;

    t = 0;
    dt = 0.005;
    s = t/tau;
    s_dot = 1/tau;
    s_ddot = 0;

    n_dof3 = 3*n_dof; % for pos, vel, accel
    
    %% --------  Init MPC  --------
    N = 10;%10; %200;    
%     dt_ = dt * (1:N).^0.9;
    dt_ = 0.02*ones(1,N); %dt;
    
    N_kernels = gmp.numOfKernels();
    
    pos_slack = 0;
    vel_slack = 0;
    accel_slack = 0;
    n_slack = pos_slack + vel_slack + accel_slack;
    
    Aineq_slack = [];
    Q_slack = [];
    if (pos_slack)
        Q_slack = blkdiag(Q_slack, 1000); 
        Aineq_slack = [Aineq_slack [1; 0; 0]];
    end
    if (vel_slack)
        Q_slack = blkdiag(Q_slack, 100);
        Aineq_slack = [Aineq_slack [0; 1; 0]];
    end
    if (accel_slack)
        Q_slack = blkdiag(Q_slack, 10);
        Aineq_slack = [Aineq_slack [0; 0; 1]];
    end
    Q_slack = sparse(Q_slack);
    Aineq_slack = sparse( repmat(Aineq_slack, n_dof, 1) );
    
    n = n_dof * N_kernels + n_slack;
    
    % state for control horizon
    % Z = [x(1), x(2), ... x(N), u(0), u(1), ... u(N-1)]
    
    % State tracking gains: (x(i) - xd(i))'*Qi*(x(i) - xd(i))
    Qi = blkdiag( opt_pos*speye(n_dof,n_dof) , opt_vel*10*speye(n_dof,n_dof));
    QN = blkdiag( 100*speye(n_dof,n_dof) , 1*speye(n_dof,n_dof));

    z_min = [pos_lim(:,1); vel_lim(:,1); accel_lim(:,1)];
    z_max = [pos_lim(:,2); vel_lim(:,2); accel_lim(:,2)];
    
    w = gmp.W';
    w = w(:); % initial guess for weights
    Z0 = [w; zeros(n_slack,1)];
    
    gmp.setTruncatedKernels(true,1e-8);

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
    Phi_f = sparse([kron(eye(n_dof),phi_f'); kron(eye(n_dof),phi_f_dot'); kron(eye(n_dof),phi_f_ddot')]);
    x_final = [yg; zeros(n_dof,1); zeros(n_dof,1)];

    %% --------  Simulation loop  --------
    while (true)
        
        %% --------  Stopping criteria  --------
        if (s > 1.0), break; end
        
        %t/tau
        text_prog.update(100*t/tau);
%         fprintf('progress: %.1f\n',t/tau);
        
        if (s >= 1)
            s = 1;
            s_dot = 0;
            s_ddot = 0;
        end

        %% -------  calc problem matrices  --------
        
%         % restore bounds in case they were changed in previous iter for si>=1
%         Z_min = repmat(z_min, N,1);
%         Z_max = repmat(z_max, N,1);
        
        H = sparse(n,n);
        q = zeros(n,1);
        Aineq = sparse(N*n_dof3 + n_slack, n);
        
        % DMP phase variable
        si = s;
        si_dot = s_dot;
        si_ddot = s_ddot;
        
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

            if (i==N), Qi_ = QN;
            else, Qi_ = Qi;
            end
            
            Psi = [kron(speye(n_dof),phi'); kron(speye(n_dof),phi_dot')];
            xd_i = [yd_i; dyd_i];

            H = H + blkdiag(Psi'*Qi_*Psi, Q_slack);
            q = q - [Psi'*Qi_*xd_i; zeros(n_slack,1)];
            
            Aineq_i = [kron(speye(n_dof),phi'); kron(speye(n_dof),phi_dot'); kron(speye(n_dof),phi_ddot')];
            Aineq((i-1)*n_dof3+1 : i*n_dof3, :) = [Aineq_i, Aineq_slack];
            
            si = si + si_dot*dt_(i);
            si_dot = si_dot + si_ddot*dt_(i);
            % si_ddot = ... (if it changes too)
        end
                      
        H = 1e-6*speye(n) + (H+H')/2; % to account for numerical errors
              
%         size(H,1)*size(H,2)
%         nnz(H)
%         
%         size(Aineq,1)*size(Aineq,2)
%         nnz(Aineq)
%         
%         Aineq
%         H
%         
%         pause
        
        Z_min = [repmat(z_min, N,1); zeros(n_slack,1)];
        Z_max = [repmat(z_max, N,1); inf(n_slack,1)];

        % This is required when horizon N shrinks towards the end
        Aineq = sparse(Aineq(1:N*n_dof3 + n_slack,:));
        Aineq(end-n_slack+1:end, end-n_slack+1:end) = speye(n_slack);
        
        Psi0 = [kron(speye(n_dof),phi0'); kron(speye(n_dof),phi0_dot'); kron(speye(n_dof),phi0_ddot')];
        Aeq = [Psi0, sparse(n_dof3, n_slack)];
        beq = [y0_; y0_dot; y0_ddot];
        
        if (si >= 1)
            Aeq = [Aeq; [Phi_f, sparse(n_dof3, n_slack)] ];
            beq = [beq; x_final];
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
            
            [Z, ~, ex_flag, opt_output] = quadprog(H,q, A,b, Aeq,beq, [],[], Z0, solver_opt);
            
            if (ex_flag == 1 || ex_flag == 2)
                % success
            else
                warning(opt_output.message);
                text_prog.printInNewLine();
                if (ex_flag < 0), break; end
            end
            
            Z0 = Z;
            w = Z(1:end-n_slack);
            slack_var = Z(end-n_slack+1:end);
            
            if (norm(slack_var) > 1e-3)
                warning('Slack activated!');
                text_prog.printInNewLine();
            end
            
            W = reshape(w, N_kernels, n_dof)';

        end
        
        %% --------  Simulate dynamics  --------
        
        yd = W*gmp.regressVec(s);
        yd_dot = W*gmp.regressVecDot(s, s_dot);
        yd_ddot = W*gmp.regressVecDDot(s, s_dot, s_ddot);

%         y_ddot = -K*(y-yd) -D*(y_dot-yd_dot) + yd_ddot;
        
        y = yd;
        y_dot = yd_dot;
        y_ddot = yd_ddot;
        
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
        y = y + y_dot*dt;
        y_dot = y_dot + y_ddot*dt;
    
    end

%     Time = [Time tau];
%     P_data = [P_data W*gmp.regressVec(1)];
%     dP_data = [dP_data W*gmp.regressVecDot(1, 0)];
%     ddP_data = [ddP_data W*gmp.regressVecDDot(1, 0, 0)];
    
    text_prog.update(100);
    fprintf('\n');
    fprintf('===> online-GMP-weights optimization finished! Elaps time: %f ms\n',toc()*1000);

end
