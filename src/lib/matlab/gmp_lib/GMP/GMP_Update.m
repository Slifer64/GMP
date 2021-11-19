%% N-DoF GMP Update class
%

classdef GMP_Update < matlab.mixin.Copyable

    methods (Access = public)
        
        %% GMP constructor.
        %  @param[in] gmp: n_DoF dmp.
        function this = GMP_Update(gmp)
                
%             if ( ~ strcmpi(class(gmp),'GMP') )
%                 error(['[GMP_Update::GMP_Update]: argument of type ''GMP'' was expected, but a type ''' class(gmp) ''' was provided instead']);
%             end
            
            this.gmp = gmp;
            
            this.initSigmaw();
            this.rv = 1;
            
            this.enableSigmawUpdate(false);

        end

        function initSigmaw(this)
            
            N_kernels = this.gmp.numOfKernels();
%             S = zeros(N_kernels,N_kernels);
%             for i=1:N_kernels
%                 for j=i+1:N_kernels
%                     S(i,j) = exp(-0.2 * abs(i-j));
%                     %S(i,j) = 1 - 1*(abs(i-j)/n)^3;
%                 end
%             end
%             S = S + S' + eye(N_kernels,N_kernels); 
%             this.Sigma_w = S;
            this.Sigma_w = eye(N_kernels, N_kernels);
            
        end
        
        function initSigmaWfromMsr(this, x_data)
              
%             n_data = length(x_data);
%             H = zeros(this.gmp.numOfKernels(), n_data);
%             H_dot = zeros(this.gmp.numOfKernels(), n_data);
%             for j=1:n_data
%                 H(:,j) = this.gmp.regressVec(x_data(j));
%                 H_dot(:,j) = this.gmp.regressVecDot(x_data(j),0.05);
%             end
%             
%             this.Sigma_w = inv(H*H' + 100*(H_dot*H_dot'));
            
            n_data = length(x_data);
            H = zeros(this.gmp.numOfKernels(), n_data);
            for j=1:n_data, H(:,j) = this.gmp.regressVec(x_data(j)); end
            
            this.Sigma_w = inv(H*H');
            
        end
        
        function enableSigmawUpdate(this, flag)
           
            this.enable_Sigma_w_update = flag;
            
        end
        
        function setSigmaW(this, Sw)
        
            this.Sigma_w = Sw;
            
        end
        
        function Sw = getSigmaW(this)
           
            Sw = this.Sigma_w;
            
        end
        
        function setMsrNoiseVar(this, rv)
           
            this.rv = rv;
            
        end
        
        %% ==================================================
        %% ===============   Online update  =================
        %% ==================================================
        
        function updatePos(this, x, y, r_n)
            
            if (nargin < 4), r_n=this.rv; end
            
            this.updateWeights(GMP_phase(x,0,0), y, GMP_UpdateType.POS, r_n);
            
        end
        
        function updateVel(this, x, x_dot, y_dot, r_n)
            
            if (nargin < 5), r_n=this.rv; end
            
            this.updateWeights(GMP_phase(x,x_dot,0), y_dot, GMP_UpdateType.VEL, r_n);
            
        end
        
        function updateAccel(this, x, x_dot, x_ddot, y_ddot, r_n)
            
            if (nargin < 6), r_n=this.rv; end
            
            this.updateWeights(GMP_phase(x,x_dot,x_ddot), y_ddot, GMP_UpdateType.ACCEL, r_n);
            
        end
        
        
        %% Updates the weights based on 'n' measurements.
        %  @param[in] s: 1 x n vector of type GMP_phase, where the j-th entry is the phase for the j-th measurement.
        %  @param[in] Z: n_dof x n matrix, where the j-th column is the j-th measurement.
        %  @param[in] type: 1 x n vector of type GMP_UpdateType, where the j-th entry has the update type for the j-th measurement.
        %  @param[in] r_n: 1 x n vector, where the j-th entry is the noise variance for the j-th measurement (optional, default=this.rv).
        function updateWeights(this, s, Z, type, r_n)
 
            if (nargin < 5), r_n = this.rv; end
                
            n = size(Z,2);
            n_ker = this.gmp.numOfKernels();
            sc = this.gmp.getScaling();
            inv_sc = this.gmp.getInvScaling();
            
            if (length(r_n) == 1), r_n = repmat(r_n(1), 1,n); end
            Rn = diag(r_n);
 
            H = zeros(n_ker, n);
            
            for j=1:n
                if (type(j) == GMP_UpdateType.POS)
                    H(:,j) = this.gmp.regressVec(s(j).x);
                    Z(:,j) = Z(:,j) - this.gmp.Y0 + sc*this.gmp.Y0d;
                elseif (type(j) == GMP_UpdateType.VEL)
                    H(:,j) = this.gmp.regressVecDot(s(j).x, s(j).x_dot);
                else % (type == GMP_UpdateType.ACCEL)
                    H(:,j) = this.gmp.regressVecDDot(s(j).x, s(j).x_dot, s(j).x_ddot);
                end
                %Z(:,j) = Z(:,j)./sc;
            end
            Z = inv_sc*Z;
            
            C = this.Sigma_w*H;
            B = eye(n,n) / (H'*this.Sigma_w*H + Rn) * C';
            this.gmp.W = this.gmp.W + (Z - this.gmp.W*H)*B;
            
            if (this.enable_Sigma_w_update)
                this.Sigma_w = this.Sigma_w - C*B;
            end
            
        end

        function plotWeightsCovariance(this, ax)
            
            if (nargin < 2)
                figure;
                ax = axes();
            end
            imshow(this.Sigma_w, [min(this.Sigma_w(:)) max(this.Sigma_w(:))], 'InitialMagnification', 2000, 'Parent',ax);
            
        end
        
    end
    
    properties (Access = protected)

        gmp % n_DoF GMP
        
        Sigma_w % weights covariance for trajectory update
        rv % noise variance for trajectory update
        enable_Sigma_w_update

    end
    
end
