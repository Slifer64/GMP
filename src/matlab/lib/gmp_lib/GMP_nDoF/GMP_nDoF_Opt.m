%% N-DoF GMP Optimization class
%

classdef GMP_nDoF_Opt < matlab.mixin.Copyable
    
    methods (Access = public)
        
        %% GMP constructor.
        %  @param[in] gmp: n_DoF dmp.
        function this = GMP_nDoF_Opt(gmp)
                
            this.ex_flag_map = containers.Map('KeyType','double','ValueType','char');
            this.ex_flag_map(-10) = 'Empty... Call ''constrOpt'' first.';
            this.ex_flag_map(1) = 'Function converged to the solution x.';
            this.ex_flag_map(0) = 'Number of iterations exceeded options.MaxIterations.';

            this.ex_flag_map(-2) = 'Problem is infeasible. Or, for interior-point-convex, the step size was smaller than options.StepTolerance, but constraints were not satisfied.';
            this.ex_flag_map(-3) = 'Problem is unbounded';
            this.ex_flag_map(2) = 'Step size was smaller than options.StepTolerance, constraints were satisfied.';
            this.ex_flag_map(-6) = 'Nonconvex problem detected.';
            this.ex_flag_map(-8) = 'Unable to compute a step direction.';
            
            this.exit_msg = '';
                
            this.gmp = gmp;
            
            this.tau = 8;
            
            this.setOptions(1,0,0, 1,1,1);
            
        end
        
        function setOptions(this, opt_pos, opt_vel, opt_accel, pos_obj_w, vel_obj_w, accel_obj_w)
            
            this.opt_pos = opt_pos;
            this.opt_vel = opt_vel;
            this.opt_accel = opt_accel;
            
            this.w_p = pos_obj_w;
            this.w_v = vel_obj_w;
            this.w_a = accel_obj_w;

        end       
        
        %% Trajectory optimization.
        %  @param[in] num_points: Number of discreet points to use in the objective function.
        %  @param[in] tau: Duration of motion.
        %  @param[in] pos_constr: Vector of @GMPConstr position constraints. For no constraints pass '[]'.
        %  @param[in] vel_constr: Vector of @GMPConstr velocity constraints. For no constraints pass '[]'.
        %  @param[in] accel_constr: Vector of @GMPConstr acceleration constraints. For no constraints pass '[]'.
        %  @param[out] success: true if minimum found, false otherwise.

        function success = optimize(this, num_points)

            if (nargin < 2), num_points = 200; end
            
            n_ker = this.gmp.numOfKernels();
            n_dof = this.gmp.numOfDoFs();
            
            x_dot = 1/this.tau;
            x_ddot = 0;
            
            success = true;
            this.exit_msg = '';
            
            x_data = linspace(0,1, num_points);
            
            % calculate cost function: J = 0.5w'Hw + f'w
            N = length(x_data);
            
            H = 1e-8*eye(n_ker, n_ker); % for numerical stability
            f = zeros(n_ker, n_dof);
            
            if (this.opt_pos)
                y_offset = - this.gmp.getY0() + this.gmp.getScaling()*this.gmp.getY0d();
                H1 = zeros(n_ker, n_ker);
                f1 = zeros(n_ker, n_dof);
                for i=1:N
                    phi = this.gmp.regressVec(x_data(i));
                    H1 = H1 + phi*phi';
                    yd = this.gmp.getYd(x_data(i)) + y_offset;
                    f1 = f1 - phi*yd';
                end
                H = H + this.w_p*H1;
                f = f + this.w_p*f1;
            end
  
            if (this.opt_vel)
                H2 = zeros(n_ker, n_ker);
                f2 = zeros(n_ker, 1);
                for i=1:N
                    phi = this.gmp.regressVecDot(x_data(i), x_dot);
                    H2 = H2 + phi*phi';
                    f2 = f2 - phi*this.gmp.getYdDot(x_data(i), x_dot)';
                end
                H = H + this.w_v*H2;
                f = f + this.w_v*f2;
            end
            
            if (this.opt_accel)
                H3 = zeros(n_ker, n_ker);
                f3 = zeros(n_ker, 1);
                for i=1:N
                    phi = this.gmp.regressVecDDot(x_data(i), x_dot, x_ddot);
                    H3 = H3 + phi*phi';
                    f3 = f3 - phi*this.gmp.getYdDDot(x_data(i), x_dot, x_ddot)';
                end
                H = H + this.w_a*H3;
                f = f + this.w_a*f3;
            end


            % inequality constraints
            
            A = [];
            b = [];

            if (~isempty(this.A_p))
                A = [A; -this.A_p; this.A_p];
                b = [b; -this.pos_lb; this.pos_ub];
            end
            
            if (~isempty(this.A_v))
                A = [A; -this.A_v; this.A_v];
                b = [b; -this.vel_lb; this.vel_ub];
            end
            
            if (~isempty(this.A_a))
                A = [A; -this.A_a; this.A_a];
                b = [b; -this.accel_lb; this.accel_ub];
            end
            
            % equality constraints
            
            Aeq = [];
            beq = [];
            
            if (~isempty(this.Aeq_p))
                Aeq = [Aeq; this.Aeq_p];
                beq = [beq; this.pos_eq];
            end
            
            if (~isempty(this.Aeq_v))
                Aeq = [Aeq; this.Aeq_v];
                beq = [beq; this.vel_eq];
            end
            
            if (~isempty(this.Aeq_a))
                Aeq = [Aeq; this.Aeq_a];
                beq = [beq; this.accel_eq];
            end
        
            % H = sparse(H);
            % A = sparse(A);
            
            opt = optimoptions('quadprog', 'Algorithm','interior-point-convex', 'StepTolerance',0, 'Display','off');
            
            % solve optimization problem
            % tic
            W = zeros(n_dof, n_ker);
            for i=1:n_dof
                
                bi = [];
                beq_i = [];
                if (~isempty(b)), bi = b(:,i); end
                if (~isempty(beq)), beq_i = beq(:,i); end
                
                [W(i,:), ~, ex_flag] = quadprog(H,f(:,i), A,bi, Aeq,beq_i, [],[], this.gmp.W(i,:), opt);

                if (ex_flag == 1 || ex_flag == 2)
                    %this.gmp.W(i,:) = w'; %w'/ks(i);
                else
                    success = false;
                    % break; ?
                end
                
                if (~isempty(this.exit_msg)), this.exit_msg = [this.exit_msg '\n']; end
                this.exit_msg = [this.exit_msg 'DoF-' num2str(i) ': ' this.ex_flag_map(ex_flag)];

            end
            
            if (success)
                this.gmp.W = this.gmp.getInvScaling()*W;
            end
            
            % toc
            
        end
                
        function setMotionDuration(this, tau)
           
            this.tau = tau;
            
        end
        
        function setPosConstr(this, x, lb, ub, x_eq, p_eq)

            n_ker = this.gmp.numOfKernels();
            n_dof = this.gmp.numOfDoFs();
            
            y_offset = - this.gmp.getY0() + this.gmp.getScaling()*this.gmp.getY0d();
            
            % Process inequality constraints
            if (~isempty(x))
 
                m = length(x);
                
                if (size(lb,2) ~= m), error(['[GMP_nDoF_Opt::setPosConstr]: Lower bounds must have ' num2str(m) ' columns (constraints).']); end
                if (size(ub,2) ~= m), error(['[GMP_nDoF_Opt::setPosConstr]: Upper bounds must have ' num2str(m) ' columns (constraints).']); end
                if (size(lb,1) ~= size(ub,1)), error('[GMP_nDoF_Opt::setPosConstr]: Lower and Upper bounds must have the same number of rows (DoFs)'); end
                
                this.A_p = zeros(m, n_ker);
                this.pos_lb = zeros(m, n_dof);
                this.pos_ub = zeros(m, n_dof);

                for i=1:m
                    this.pos_lb(i,:) = ( lb(:,i) + y_offset )';
                    this.pos_ub(i,:) = ( ub(:,i) + y_offset )';
                    this.A_p(i,:) = this.gmp.regressVec(x(i))';
                end
            
            end
            
            % Process equality constraints
            if (~isempty(x_eq))
                
                m = length(x_eq);
                if (size(p_eq,2) ~= m), error(['[GMP_nDoF_Opt::setPosConstr]: Equality values matrix must have ' num2str(m) ' columns (constraints).']); end
                
                this.Aeq_p = zeros(m, n_ker);
                this.pos_eq = zeros(m, n_dof);

                for i=1:m
                    this.pos_eq(i,:) = ( p_eq(:,i) + y_offset )';
                    this.Aeq_p(i,:) = this.gmp.regressVec(x_eq(i))';
                end
            
            end
            
        end
        
        function setVelConstr(this, x, lb, ub, x_eq, v_eq)
            
            n_ker = this.gmp.numOfKernels();
            x_dot = 1 / this.tau;

            % Process inequality constraints
            if (~isempty(x))
 
                m = length(x);
                if (size(lb,2) ~= m), error(['[GMP_nDoF_Opt::setVelConstr]: Lower bounds must have ' num2str(m) ' columns (constraints).']); end
                if (size(ub,2) ~= m), error(['[GMP_nDoF_Opt::setVelConstr]: Upper bounds must have ' num2str(m) ' columns (constraints).']); end
                if (size(lb,1) ~= size(ub,1)), error('[GMP_nDoF_Opt::setVelConstr]: Lower and Upper bounds must have the same number of rows (DoFs)'); end
                
                this.vel_lb = lb';
                this.vel_ub = ub';
                this.A_v = zeros(m, n_ker);
                for i=1:m, this.A_v(i,:) = this.gmp.regressVecDot(x(i), x_dot)'; end
            
            end
            
            % Process equality constraints
            if (~isempty(x_eq))
                
                m = length(x_eq);
                if (size(v_eq,2) ~= m), error(['[GMP_nDoF_Opt::setVelConstr]: Equality values matrix must have ' num2str(m) ' columns (constraints).']); end

                this.vel_eq = v_eq';
                this.Aeq_v = zeros(m, n_ker);
                for i=1:m, this.Aeq_v(i,:) = this.gmp.regressVecDot(x_eq(i),x_dot)'; end
            
            end
            
        end
        
        function setAccelConstr(this, x, lb, ub, x_eq, a_eq)
            
            n_ker = this.gmp.numOfKernels();
            x_dot = 1 / this.tau;
            x_ddot = 0;

            % Process inequality constraints
            if (~isempty(x))
 
                m = length(x);
                if (size(lb,2) ~= m), error(['[GMP_nDoF_Opt::setAccelConstr]: Lower bounds must have ' num2str(m) ' columns (constraints).']); end
                if (size(ub,2) ~= m), error(['[GMP_nDoF_Opt::setAccelConstr]: Upper bounds must have ' num2str(m) ' columns (constraints).']); end
                if (size(lb,1) ~= size(ub,1)), error('[GMP_nDoF_Opt::setAccelConstr]: Lower and Upper bounds must have the same number of rows (DoFs)'); end
                
                this.accel_lb = lb';
                this.accel_ub = ub';
                this.A_a = zeros(m, n_ker);
                for i=1:m, this.A_a(i,:) = this.gmp.regressVecDDot(x(i), x_dot, x_ddot)'; end
            
            end
            
            % Process equality constraints
            if (~isempty(x_eq))
                
                m = length(x_eq);
                if (size(a_eq,2) ~= m), error(['[GMP_nDoF_Opt::setAccelConstr]: Equality values matrix must have ' num2str(m) ' columns (constraints).']); end

                this.accel_eq = a_eq';
                this.Aeq_a = zeros(m, n_ker);
                for i=1:m, this.Aeq_a(i,:) = this.gmp.regressVecDDot(x_eq(i), x_dot, x_ddot)'; end
            
            end
            
        end
        
        function clearPosConstr(this)
                
            this.A_p = [];
            this.pos_lb = [];
            this.pos_ub = [];
            
            this.Aeq_p = [];
            this.pos_eq = [];

        end
        
        function clearVelConstr(this)
                
            this.A_v = [];
            this.vel_lb = [];
            this.vel_ub = [];
            
            this.Aeq_v = [];
            this.vel_eq = [];
            
        end
        
        function clearAccelConstr(this)
            
            this.A_a = [];
            this.accel_lb = [];
            this.accel_ub = [];
            
            this.Aeq_a = [];
            this.accel_eq = [];
            
        end
        
        function clearConstr(this)
    
            this.clearPosConstr();
            this.clearVelConstr();
            this.clearAccelConstr();
            
        end
        
        function msg = getExitMsg(this)
           
            msg = this.exit_msg;
            
        end

    end
    
    properties (Access = private)
        
        ex_flag_map % maps the exit flag value to the msg describing the exit flag
        exit_msg
        
        gmp % n_DoF GMP
        
        opt_pos % flag indicating to optimize position
        opt_vel % flag indicating to optimize velocity
        opt_accel % flag indicating to optimize acceleration
        
        w_p % position objective weight
        w_v % velocity objective weight
        w_a % acceleration objective weight
        
        tau % motion duration
        
        % Constraints
        A_p
        pos_lb
        pos_ub
        
        Aeq_p
        pos_eq
        
        A_v
        vel_lb
        vel_ub
        
        Aeq_v
        vel_eq
        
        A_a
        accel_lb
        accel_ub
        
        Aeq_a
        accel_eq

    end
    
end
