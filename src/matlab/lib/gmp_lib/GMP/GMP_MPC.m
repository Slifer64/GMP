classdef GMP_MPC < handle
   
    methods (Access = public)
       
        function this = GMP_MPC(gmp, N_horizon, pred_time_step, N_kernels, kernel_std_scaling, slack_flags, slack_gains)

            n_dof = gmp.numOfDoFs();
            
            this.gmp_ref = gmp;
            
            this.gmp_mpc = GMP(n_dof, N_kernels, kernel_std_scaling);
            this.gmp_mpc.setTruncatedKernels(true,1e-8);
            
            this.N = N_horizon;
            this.dt_ = pred_time_step*ones(1,this.N);
  
            opt_pos = 1;
            opt_vel = 0;
            this.pos_slack = slack_flags(1);
            this.vel_slack = slack_flags(2);
            this.accel_slack = slack_flags(3);
            this.n_slack = this.pos_slack + this.vel_slack + this.accel_slack;

            this.Aineq_slack = [];
            this.Q_slack = [];
            if (this.pos_slack)
                this.Q_slack = blkdiag(this.Q_slack, slack_gains(1)); 
                this.Aineq_slack = [this.Aineq_slack [-ones(n_dof,1); zeros(2*n_dof,1)] ];
            end
            if (this.vel_slack)
                this.Q_slack = blkdiag(this.Q_slack, slack_gains(2));
                this.Aineq_slack = [this.Aineq_slack [zeros(n_dof,1); -ones(n_dof,1); zeros(n_dof,1)] ];
            end
            if (this.accel_slack)
                this.Q_slack = blkdiag(this.Q_slack, slack_gains(3));
                this.Aineq_slack = [this.Aineq_slack [zeros(2*n_dof,1); -ones(n_dof,1)] ];
            end
            this.Q_slack = sparse(this.Q_slack);
            this.Aineq_slack = sparse(this.Aineq_slack);

            % State tracking gains: (x(i) - xd(i))'*Qi*(x(i) - xd(i))
            this.Qi = blkdiag( opt_pos*speye(n_dof,n_dof) , opt_vel*10*speye(n_dof,n_dof));
            
            % if I set it to Qi or zero the velocity part, sometimes it fails to find feasible solution 
            this.QN = blkdiag( 100*speye(n_dof,n_dof) , 1*speye(n_dof,n_dof) ); %this.Qi;

            n_dof3 = 3*n_dof;
            this.Z0 = zeros(n_dof*N_kernels + this.n_slack,1);
            this.Z0_dual_ineq = zeros(this.N*n_dof3 + this.n_slack, 1);
            this.Z0_dual_eq = zeros(n_dof3, 1);
            
            this.setPosLimits(-inf(n_dof,1), inf(n_dof,1));
            this.setVelLimits(-inf(n_dof,1), inf(n_dof,1));
            this.setAccelLimits(-inf(n_dof,1), inf(n_dof,1));
            
            this.setPosSlackLimit(inf);
            this.setVelSlackLimit(inf);
            this.setAccelSlackLimit(inf);

        end
        
        function setInitialState(this, y0, y0_dot, y0_ddot, s, s_dot, s_ddot)

            n_dof = this.gmp_mpc.numOfDoFs();
            phi0 = this.gmp_mpc.regressVec(s);
            phi0_dot = this.gmp_mpc.regressVecDot(s, s_dot);
            phi0_ddot = this.gmp_mpc.regressVecDDot(s, s_dot, s_ddot);
            this.Phi0 = sparse([kron(eye(n_dof),phi0'); kron(eye(n_dof),phi0_dot'); kron(eye(n_dof),phi0_ddot')]);
            this.x0 = [y0; y0_dot; y0_ddot];
            
        end
        
        function setFinalState(this, yf, yf_dot, yf_ddot, s, s_dot, s_ddot)

            n_dof = this.gmp_mpc.numOfDoFs();
            this.s_f = s;
            phi_f = this.gmp_mpc.regressVec(s);
            phi_f_dot = this.gmp_mpc.regressVecDot(s, s_dot);
            phi_f_ddot = this.gmp_mpc.regressVecDDot(s, s_dot, s_ddot);
            this.Phi_f = sparse([kron(eye(n_dof),phi_f'); kron(eye(n_dof),phi_f_dot'); kron(eye(n_dof),phi_f_ddot')]);
            this.x_f = [yf; yf_dot; yf_ddot];
            
        end

        function setPosLimits(this, lb, ub)
            
            this.pos_lb = lb;
            this.pos_ub = ub;
            
        end
        
        function setVelLimits(this, lb, ub)
            
            this.vel_lb = lb;
            this.vel_ub = ub;
            
        end
        
        function setAccelLimits(this, lb, ub)
            
            this.accel_lb = lb;
            this.accel_ub = ub;
            
        end
        
        function setPosSlackLimit(this, s_lim)
            
            this.pos_slack_lim = s_lim;
        end
        
        function setVelSlackLimit(this, s_lim)
            
            this.vel_slack_lim = s_lim;
        end
        
        function setAccelSlackLimit(this, s_lim)
            
            this.accel_slack_lim = s_lim;
        end
  
        function [exit_flag, y, y_dot, y_ddot, slack_var] = solve(this, s, s_dot, s_ddot)
            
            y = [];
            y_dot = []; 
            y_ddot = [];
            slack_var = [];
            
            N_kernels = this.gmp_mpc.numOfKernels();
            n_dof = this.gmp_ref.numOfDoFs();
            n_dof3 = 3*n_dof; % for pos, vel, accel
            
            n = n_dof * N_kernels + this.n_slack;
            
            H = 1e-6*speye(n); % for numerical stability
            q = zeros(n,1);
            % add slacks to bounds. I.e. one could optionaly define with slacks
            % bounds the maximum allowable error
            Aineq = sparse(this.N*n_dof3 + this.n_slack, n);
            Aineq(end-this.n_slack+1:end,end-this.n_slack+1:end) = speye(this.n_slack);
            
            z_min = [this.pos_lb; this.vel_lb; this.accel_lb];
            z_max = [this.pos_ub; this.vel_ub; this.accel_ub];
            
            slack_lim = [];
            if (this.pos_slack), slack_lim = [slack_lim; this.pos_slack_lim];  end
            if (this.vel_slack), slack_lim = [slack_lim; this.vel_slack_lim];  end
            if (this.accel_slack), slack_lim = [slack_lim; this.accel_slack_lim];  end
            
            Z_min = [repmat(z_min, this.N,1); -slack_lim];
            Z_max = [repmat(z_max, this.N,1); slack_lim];

            % DMP phase variable
            si = s;
            si_dot = s_dot;
            si_ddot = s_ddot;

            for i=1:this.N

                yd_i = this.gmp_ref.getYd(si);
                dyd_i = this.gmp_ref.getYdDot(si, si_dot);

                phi = this.gmp_mpc.regressVec(si);
                phi_dot = this.gmp_mpc.regressVecDot(si, si_dot);
                phi_ddot = this.gmp_mpc.regressVecDDot(si, si_dot, si_ddot);

                if (i==this.N), Qi_ = this.QN;
                else, Qi_ = this.Qi;
                end

                Psi = [kron(speye(n_dof),phi'); kron(speye(n_dof),phi_dot')];
                xd_i = [yd_i; dyd_i];

                H = H + blkdiag(Psi'*Qi_*Psi, this.Q_slack);
                q = q - [Psi'*Qi_*xd_i; zeros(this.n_slack,1)];

                Aineq_i = [kron(speye(n_dof),phi'); kron(speye(n_dof),phi_dot'); kron(speye(n_dof),phi_ddot')];
                Aineq((i-1)*n_dof3+1 : i*n_dof3, :) = [Aineq_i, this.Aineq_slack];

                si = si + si_dot*this.dt_(i);
                si_dot = si_dot + si_ddot*this.dt_(i);
                % si_ddot = ... (if it changes too)
            end

            H = (H+H')/2; % to account for numerical errors

            Aeq = [this.Phi0, sparse(n_dof3, this.n_slack)];
            beq = [this.x0];

            if (si >= this.s_f)
                Aeq = [Aeq; [this.Phi_f, sparse(n_dof3, this.n_slack)] ];
                beq = [beq; this.x_f];

                n_eq_plus = size(Aeq,1) - length(this.Z0_dual_eq);
                if (n_eq_plus>0), this.Z0_dual_eq = [this.Z0_dual_eq; zeros(n_eq_plus, 1)]; end
            end

            %% ===========  solve optimization problem  ==========

            A_osqp = [Aineq; Aeq];
            lb = [Z_min; beq];
            ub = [Z_max; beq];

            Z0_dual = [this.Z0_dual_ineq; this.Z0_dual_eq];

            % Create an OSQP object
            osqp_solver = osqp;
            osqp_solver.setup(H, q, A_osqp, lb, ub, 'warm_start',false, 'verbose',false, 'eps_abs',1e-4, 'eps_rel',1e-5);%, 'max_iter',20000);
            osqp_solver.warm_start('x', this.Z0, 'y',Z0_dual);

            res = osqp_solver.solve();

            this.exit_msg = ''; % clear exit message
            exit_flag = 0;
            if ( res.info.status_val ~= 1)
                %res.info
                this.exit_msg = res.info.status;
                exit_flag = 1;
                if (res.info.status_val == -3 || res.info.status_val == -4 || res.info.status_val == -7 || res.info.status_val == -10)
                    exit_flag = -1;
                    return; 
                end
            end

            Z = res.x;
            this.Z0 = Z;
            Z0_dual = res.y;
            n_ineq = size(Aineq,1);
            this.Z0_dual_ineq = Z0_dual(1:n_ineq);
            this.Z0_dual_eq = Z0_dual(n_ineq+1:end);
            
            w = Z(1:end-this.n_slack);
            slack_var = Z(end-this.n_slack+1:end);
            W = reshape(w, N_kernels, n_dof)';
     
            %% --------  Generate output  --------

            y = W*this.gmp_mpc.regressVec(s);
            y_dot = W*this.gmp_mpc.regressVecDot(s, s_dot);
            y_ddot = W*this.gmp_mpc.regressVecDDot(s, s_dot, s_ddot);

            this.setInitialState(y, y_dot, y_ddot, s, s_dot, s_ddot);
            
        end
        
        function msg = getExitMsg(this)
           
            msg = this.exit_msg;
            
        end
        
    end
    
    properties (Access = protected)
        
        exit_msg % std::string
        
        gmp_ref % const GMP*

        gmp_mpc % std::shared_ptr<GMP>

        N % unsigned
        dt_ % arma::rowvec(N)
        
        Aineq_slack % arma::mat(3*n_dof, n_slack)
        Q_slack % arma::mat(n_slack, n_slack)
        
        Qi % arma::mat(n_dof*N_kernels, n_dof*N_kernels)
        QN % arma::mat(n_dof*N_kernels, n_dof*N_kernels)
        
        Z0 
        Z0_dual_ineq
        Z0_dual_eq
        
        pos_lb % arma::vec(n_dof)
        pos_ub % arma::vec(n_dof)
        
        vel_lb % arma::vec(n_dof)
        vel_ub % arma::vec(n_dof)
        
        accel_lb % arma::vec(n_dof)
        accel_ub % arma::vec(n_dof)
        
        n_slack % unsigned
        
        pos_slack % bool
        vel_slack % bool
        accel_slack % bool
        
        pos_slack_lim % double
        vel_slack_lim % double
        accel_slack_lim % double
        
        Phi0 % arma::mat(3*n_dof, n_dof*N_kernels)
        x0 % arma::vec(3*n_dof)
        
        s_f % double
        Phi_f % arma::mat(3*n_dof, n_dof*N_kernels)
        x_f % arma::vec(3*n_dof)
        
    end
    
    
end