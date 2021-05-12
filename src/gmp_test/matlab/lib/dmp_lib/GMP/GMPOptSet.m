%% GMP optimization options
%

classdef GMPOptSet < matlab.mixin.Copyable
    
    properties (Access = public, Constant)

        
    end
    
    methods (Access = public)
        
        %% Constructor.
        %  @param[in] opt_pos: flag indicating to optimize position
        %  @param[in] opt_vel: flag indicating to optimize velocity
        %  @param[in] opt_accel: flag indicating to optimize acceleration
        %  @param[in] pos_obj_w: position objective weight
        %  @param[in] vel_obj_w: velocity objective weight
        %  @param[in] accel_obj_w: acceleration objective weight
        function this = GMPOptSet(opt_pos, opt_vel, opt_accel, pos_obj_w, vel_obj_w, accel_obj_w)
                
            if (nargin < 6), return; end
            
            this.opt_pos = opt_pos;
            this.opt_vel = opt_vel;
            this.opt_accel = opt_accel;
            
            this.w_p = pos_obj_w;
            this.w_v = vel_obj_w;
            this.w_a = accel_obj_w;

        end

    end
    
    properties (Access = public)
     
        opt_pos % flag indicating to optimize position
        opt_vel % flag indicating to optimize velocity
        opt_accel % flag indicating to optimize acceleration
        
        w_p % position objective weight
        w_v % velocity objective weight
        w_a % acceleration objective weight
        
    end
    
end
