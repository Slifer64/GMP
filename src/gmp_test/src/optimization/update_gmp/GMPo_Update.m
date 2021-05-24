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
        
        % function setSigmaW(this, Sw)
        
        % function Sw = getSigmaW(this)
        
        % function setMsrNoiseVar(this, rv)
        
        
        %% ==================================================
        %% ===============   Online update  =================
        %% ==================================================
        
        % function updatePos(this, x, y, r_n)
        
        % function updateVel(this, x, x_dot, y_dot, r_n)
        
        % function updateAccel(this, x, x_dot, x_ddot, y_ddot, r_n)
        
        function updateQuat(this, x, Q, r_n)
            
            if (nargin < 4), r_n=this.rv; end
            
            y = GMPo.quat2q(Q, this.gmp.getQ0());
            
            this.updateWeights(GMP_phase(x,0,0), y, GMPo_UpdateType.POS, r_n);
            
        end
        
        function updateRotVel(this, x, x_dot, rotVel, Q, r_n)
            
            if (nargin < 6), r_n=this.rv; end
            
            Q1 = GMPo.quat2q(Q, this.gmp.getQ0());
            y_dot = gmp_.rotVel_to_qLogDot(rotVel, Q1);
            
            this.updateWeights(GMP_phase(x,x_dot,0), y_dot, GMPo_UpdateType.VEL, r_n);
            
        end
        
        function updateRotAccel(this, x, x_dot, x_ddot, rotAccel, rotVel, Q, r_n)
            
            if (nargin < 6), r_n=this.rv; end
            
            Q1 = GMPo.quat2q(Q, this.gmp.getQ0());
            y_ddot = rotAccel_to_qLogDDot(rotAccel, rotVel, Q1);
            
            this.updateWeights(GMP_phase(x,x_dot,x_ddot), y_ddot, GMPo_UpdateType.ACCEL, r_n);
            
        end

        % function updateWeights(this, s, Z, type, r_n)
        
        % function plotWeightsCovariance(this, ax)
    end
    
end
