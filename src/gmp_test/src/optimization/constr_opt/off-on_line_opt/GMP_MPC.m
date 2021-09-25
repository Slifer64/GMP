classdef GMP_MPC
   
    
    methods (Access = public)
       
        
        function this = GMP_MPC(gmp, N_kernels, kernel_std_scaling)
            
            if (nargin < 2), N_kernels = gmp.numOfKernels(); end
            if (nargin < 3), kernel_std_scaling = gmp.numOfKernels(); end
            
            this.gmp_ = gmp;
           
            this.n_dof = length(gmp.numOfDoFs());

            N_kernels = 30; %gmp0.numOfKernels();
            s_data = 0:0.01:1;
            Yd_data = zeros(n_dof, length(s_data));
            for j=1:length(s_data), Yd_data(:,j)=gmp.getYd(s_data(j)); end
            gmp = GMP(n_dof, N_kernels, 1.5);
            gmp.train('LS', s_data, Yd_data);

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
            N = 10; %15;    
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
                Q_slack = blkdiag(Q_slack, 1000000); 
                Aineq_slack = [Aineq_slack [-ones(n_dof,1); zeros(2*n_dof,1)] ];
            end
            if (vel_slack)
                Q_slack = blkdiag(Q_slack, 100);
                Aineq_slack = [Aineq_slack [zeros(n_dof,1); -ones(n_dof,1); zeros(n_dof,1)] ];
            end
            if (accel_slack)
                Q_slack = blkdiag(Q_slack, 20);
                Aineq_slack = [Aineq_slack [zeros(2*n_dof,1); -ones(n_dof,1)] ];
            end
            Q_slack = sparse(Q_slack);
            Aineq_slack = sparse(Aineq_slack);

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
            n_ineq = N*n_dof3 + n_slack;
            n_eq = n_dof3;
            Z0_dual_ineq = zeros(n_ineq, 1);
            Z0_dual_eq = zeros(n_eq, 1);

            gmp.setTruncatedKernels(true,1e-8);

            %% --------  Init solver  --------
            if (qp_solver_type == 0)
                solver_opt = optimoptions('quadprog', 'Algorithm','interior-point-convex', 'LinearSolver','sparse', 'StepTolerance',0, 'Display','off', 'MaxIterations',2000);
            end

            text_prog = ProgressText(40);
            text_prog.init();

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
        
        end
        
        
    end
    
    properties (Access = protected)
        
        
    end
    
    
end