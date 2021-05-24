%% N-DoF GMP Update class
%

classdef GMPo_Update < GMP_Update

    methods (Access = public)
        
        %% GMP constructor.
        %  @param[in] gmp: n_DoF dmp.
        function this = GMPo_Update(gmp)
                
            this@GMP_Update(gmp);

        end

        % function initSigmaw(this)
        
        % function initSigmaWfromMsr(this, x_data)
        
        % function enableSigmawUpdate(this, flag)
        
        %function setSigmaW(this, Sw)
        
        %function Sw = getSigmaW(this)
        
        % function setMsrNoiseVar(this, rv)
        
        
        %% ==================================================
        %% ===============   Online update  =================
        %% ==================================================
        
        function updateQuat(this, x, y, r_n)
            
            if (nargin < 4), r_n=this.rv; end
            
            y = GMPo.
            
            this.updateWeights(GMP_phase(x,0,0), y, GMPo_UpdateType.POS, r_n);
            
        end
        
        function updateRotVel(this, x, x_dot, y_dot, r_n)
            
            if (nargin < 5), r_n=this.rv; end
            
            this.updateWeights(GMP_phase(x,x_dot,0), y_dot, GMPo_UpdateType.VEL, r_n);
            
        end
        
        function updateRotAccel(this, x, x_dot, x_ddot, y_ddot, r_n)
            
            if (nargin < 6), r_n=this.rv; end
            
            this.updateWeights(GMP_phase(x,x_dot,x_ddot), y_ddot, GMPo_UpdateType.ACCEL, r_n);
            
        end
        
        
        %% Updates the weights based on 'n' measurements.
        %  @param[in] s: 1 x n vector of type GMP_phase, where the j-th entry is the phase for the j-th measurement.
        %  @param[in] Z: n_dof x n matrix, where the j-th column is the j-th measurement.
        %  @param[in] type: 1 x n vector of type GMPo_UpdateType, where the j-th entry has the update type for the j-th measurement.
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
                if (type(j) == GMPo_UpdateType.POS)
                    H(:,j) = this.gmp.regressVec(s(j).x);
                    Z(:,j) = Z(:,j) - this.gmp.Y0 + sc*this.gmp.Y0d;
                elseif (type(j) == GMPo_UpdateType.VEL)
                    H(:,j) = this.gmp.regressVecDot(s(j).x, s(j).x_dot);
                else % (type == GMPo_UpdateType.ACCEL)
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
    
    properties (Access = private)

        gmp % n_DoF GMP
        
        Sigma_w % weights covariance for trajectory update
        rv % noise variance for trajectory update
        enable_Sigma_w_update

    end
    
end
