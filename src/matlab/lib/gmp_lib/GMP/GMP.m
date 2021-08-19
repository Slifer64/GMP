%% N-DoF GMP class
%  Generalized movement primitive.
%

classdef GMP < GMP_regressor
    
    methods (Access = public)
        
        %% GMP constructor.
        %  @param[in] n_dofs: number of degrees of freedom.
        %  @param[in] N_kernels: the number of kernels
        %  @param[in] kern_std_scale: Scaling for std of kernels (optional, default=1).
        function this = GMP(n_dofs, N_kernels, kern_std_scale)
                
            if (nargin < 1), n_dofs = 1; end
            if (nargin < 2), N_kernels = 2; end
            if (nargin < 3), kern_std_scale = 1.0; end
            
            this@GMP_regressor(N_kernels, kern_std_scale);

            this.K = 100 * ones(n_dofs,1);
            this.D = 30 * ones(n_dofs,1);
            
            this.W = zeros(n_dofs, N_kernels);
            
            this.Y0d = zeros(n_dofs,1);
            this.Ygd = ones(n_dofs,1);
            this.Y0 = this.Y0d;
            this.Yg = this.Ygd;
            
            this.setScaleMethod(TrajScale_Prop(n_dofs));

            this.y_dot = zeros(n_dofs,1);
            this.z_dot = zeros(n_dofs,1);
            
        end

        %% Returns the number of DoFs.
        %  return: number of DoFs.
        function n = numOfDoFs(this)
           
            n = size(this.W,1);
            
        end
        
        %% Returns the number of kernels.
        %  @return: number of kernels.
        function n_ker = numOfKernels(this)
            
            n_ker = length(this.c);
            
        end

        %% Sets the initial position.
        %  @param[in] y0: initial position.
        function setY0(this, y0)
            
            this.Y0 = y0; 
            this.traj_sc.setNewStartFinalPos(this.Y0, this.Yg);
        
        end
        
        function Y0 = getY0(this), Y0 = this.Y0; end
        
        function y0d = getY0d(this), y0d = this.Y0d; end
        
        
        %% Set goal position.
        %  @param[in] g: goal position.
        function setGoal(this, g)
            
            this.Yg = g;
            this.traj_sc.setNewStartFinalPos(this.Y0, this.Yg);
        
        end

        function Yg = getGoal(this), Yg = this.Yg; end
        
        function sc = getScaling(this)

            sc = this.traj_sc.getScaling();
            
        end
        
        function inv_sc = getInvScaling(this)
            
            inv_sc = this.traj_sc.getInvScaling();
            
        end

        %% Set the trajectory spatial scaling method.
        %  @param[in] traj_scale_obj: pointer to an object of type @TrajScale.
        function setScaleMethod(this, traj_scale_obj)
           
            if (this.numOfDoFs() ~= traj_scale_obj.numOfDoFs())
                error('[GMP::setScaleMethod]: Incompatible number of DoFs...');
            end
            
            this.traj_sc = traj_scale_obj;
            this.traj_sc.setNominalStartFinalPos(this.Y0d, this.Ygd);
            this.traj_sc.setNewStartFinalPos(this.Y0, this.Yg);

        end
        
        %% Trains the GMP.
        %  @param[in] train_method: the training method to use, as a string ('LWR', 'LS').
        %  @param[in] x: Row vector with the canonical timestamps (in [0 1]) of the training data points.
        %  @param[in] yd_data: Matrix with the desired potition for each DoF in each row.
        %  @param[out] train_error: The training error expressed as the mse error.
        %  \note The timestamps-data need not be sequencial temporarily.
        function [train_err, Sw] = train(this, train_method, x, yd_data)
            
            if (~isempty(find(x>1 | x<0)))
               warning('[GMP::train]: The training timestamps are not normalized...') ;
            end
            
            n_data = length(x);
            num_ker = this.numOfKernels();
            n_dofs = this.numOfDoFs();

            H = zeros(num_ker, n_data);
            % for j=1:n_data, Psi(:,j) = this.kernelFun(x(j)); end
            for j=1:n_data, H(:,j) = this.regressVec(x(j)); end

            if (strcmpi(train_method,'LWR') == 1)
                this.W = yd_data*H' ./ repmat(sum(H,2)',n_dofs,1);
            elseif (strcmpi(train_method,'LS') == 1)
                this.W = yd_data / H;
            else
                error('[WSoG::train]: Unsupported training method...');
            end
            
            this.Y0d = this.W*H(:,1);
            this.Ygd = this.W*H(:,end);
            
            this.setY0(this.Y0d);
            this.setGoal(this.Ygd);
            
            this.traj_sc.setNominalStartFinalPos(this.Y0d, this.Ygd);
            this.traj_sc.setNewStartFinalPos(this.getY0(), this.getGoal());

            if (nargout > 0)
                train_err = zeros(n_dofs,1);
                err_data = this.W*H - yd_data;
                for i=1:n_dofs, train_err(i) = norm(err_data(i,:)); end
            end
            
            if (nargout > 1), Sw = inv(H*H'); end
            
        end

        %% Auto-retrains the GMP modifying the kernels according to the arguments.
        %  @param[in] N_kernels: the new number of kernels.
        %  @param[in] kern_std_scale: the width of the kernels.
        %  @param[in] n_points: number of uniformly autogenerated points in [0 1] for the auto-retraining.
        %  @param[in] train_method: the training method to use, as a string ('LWR', 'LS').
        function [train_err, Sw] = autoRetrain(this, N_kernels, kern_std_scale, n_points, train_method)
            
            if (nargin < 4), n_points = 500; end
            if (nargin < 5), train_method = 'LS'; end
            
            x = linspace(0,1, n_points);
            yd_data = zeros(this.numOfDoFs(), length(x));
            for j=1:length(x), yd_data(:,j) = this.W*this.regressVec(x(j)); end
            
            this.setKernels(N_kernels, kern_std_scale);
            [train_err, Sw] = this.train(train_method, x, yd_data);
        
        end

        %% Calculates the time derivatives of the GMP's states.
        %  @param[in] s: Object of type @GMP_phase.
        %  @param[in] y: 'y' state of the GMP.
        %  @param[in] z: 'z' state of the GMP.
        %  @param[in] y_c: coupling term for the dynamical equation of the 'y' state (optional, default=0).
        %  @param[in] z_c: coupling term for the dynamical equation of the 'z' state (optional, default=0).
        function update(this, s, y, z, y_c, z_c)

            n_dofs = this.numOfDoFs();
            
            if (nargin < 5), y_c = zeros(n_dofs,1); end
            if (nargin < 6), z_c = zeros(n_dofs,1); end
            
            if (length(y_c) == 1), y_c = ones(n_dofs,1)*y_c(1); end
            if (length(z_c) == 1), z_c = ones(n_dofs,1)*z_c(1); end
            
            yd = this.getYd(s.x);
            yd_dot = this.getYdDot(s.x, s.x_dot);
            yd_ddot = this.getYdDDot(s.x, s.x_dot, s.x_ddot);
  
            this.y_dot = z + y_c;
            this.z_dot = this.K.*(yd - y) + this.D.*(yd_dot - z) + yd_ddot + z_c;

        end

        
        %% Returns the 'y' state time derivative.
        %  Call after @update.
        %  @return: time derivative of 'y' state.
        function y_dot = getYdot(this)
            
            y_dot = this.y_dot;
            
        end
        
        
        %% Returns the 'z' state time derivative.
        %  Call after @update.
        %  @return: time derivative of 'z' state.
        function z_dot = getZdot(this)
            
            z_dot = this.z_dot;
            
        end
        
        
        %% Returns the GMP's acceleration.
        %  Call after @update.
        %  @param[in] yc_dot: time derivative of 'y' state coupling (optional, default=0).
        %  @return: acceleration.
        function y_ddot = getYddot(this, yc_dot)
            
            n_dofs = this.numOfDoFs();
            if (nargin < 2), yc_dot = zeros(n_dofs,1); end
            if (length(yc_dot)==1), yc_dot = ones(n_dofs,1)*yc_dot(1); end
            
            y_ddot = this.getZdot() + yc_dot;

        end
        
        
        %% Calculates the GMP's acceleration.
        %  @param[in] s: Object of type @GMP_phase.
        %  @param[in] y: 'y' state of the GMP.
        %  @param[in] y_dot: time derivative of 'y' state.
        %  @param[in] y_c: coupling term for the dynamical equation of the 'y' state (optional, default=0).
        %  @param[in] z_c: coupling term for the dynamical equation of the 'z' state (optional, default=0).
        %  @param[in] yc_dot: time derivative of the 'y' state coupling (optional, default=0).
        %  @return: acceleration.
        function y_ddot = calcYddot(this, s, y, y_dot, y_c, z_c, yc_dot)

            n_dofs = this.numOfDoFs();
            if (nargin < 5), y_c = zeros(n_dofs,1); end
            if (nargin < 6), z_c = zeros(n_dofs,1); end
            if (nargin < 7), yc_dot = zeros(n_dofs,1); end
            
            if (length(y_c)==1), y_c = ones(n_dofs,1)*y_c(1); end
            if (length(z_c)==1), z_c = ones(n_dofs,1)*z_c(1); end
            if (length(yc_dot)==1), yc_dot = ones(n_dofs,1)*yc_dot(1); end

            yd = this.getYd(s.x);
            yd_dot = this.getYdDot(s.x, s.x_dot);
            yd_ddot = this.getYdDDot(s.x, s.x_dot, s.x_ddot);
            
            z = y_dot - y_c;
            z_dot = this.K.*(yd - y) + this.D.*(yd_dot - z) + yd_ddot + z_c;
            y_ddot = (z_dot + yc_dot);

        end
        
        %% Returns the scaled desired position.
        %  @param[in] x: phase variable.
        %  @return: scaled desired position.
        function yd = getYd(this, x)

            yd = this.getScaling()*(this.W*this.regressVec(x) - this.Y0d) + this.Y0;
            
        end
        
        
        %% Returns the scaled desired velocity.
        %  @param[in] x: phase variable.
        %  @param[in] x_dot: 1st time derivative of the phase variable.
        %  @return: scaled desired velocity.
        function yd_dot = getYdDot(this, x, x_dot)
            
            yd_dot = this.getScaling()*this.W*this.regressVecDot(x,x_dot);

        end
        
        
        %% Returns the scaled desired acceleration.
        %  @param[in] x: phase variable.
        %  @param[in] x_dot: 1st time derivative of the phase variable.
        %  @param[in] x_ddot: 2nd time derivative of the phase variable.
        %  @return: scaled desired acceleration.
        function yd_ddot = getYdDDot(this, x, x_dot, x_ddot)
            
            if (nargin < 4), x_ddot = 0; end
            yd_ddot = this.getScaling()*this.W*this.regressVecDDot(x,x_dot,x_ddot);
            
        end


        %% Returns a deep copy of this object.
        %  @return: deep copy of this object.
        function cp_obj = deepCopy(this)
            
            % Make a shallow copy of all properties
            cp_obj = this.copy();
            
            % Make a deep copy of the pointers
            cp_obj.traj_sc = this.traj_sc.deepCopy();

        end
        
    end
    
    properties (Access = public)

        traj_sc % object of type @TrajScale
        
        %% weights
        W % num_DoF x num_Kernels matrix where each row contrains the weights for each DoF

        %% impedance params
        K % num_DoF x num_DoF stiffness matrix
        D % num_DoF x num_DoF damping matrix 

    end
    
    
    properties (Access = {?GMP_Update, ?gmp_})
        
        Y0 % initial position
        Yg % target position
        
        Y0d % trained initial position
        Ygd % trained target position
        
        %% output state
        y_dot % position derivative
        z_dot % scaled velocity derivative
        
    end
    
end
