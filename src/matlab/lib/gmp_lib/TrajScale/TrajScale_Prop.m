%% TrajScale_Prop class

classdef TrajScale_Prop < TrajScale
    
    methods (Access = public)
        
        %% Constructor.
        %  @param[in] n_dof: degrees of freedom.
        function this = TrajScale_Prop(n_dof)
    
            this@TrajScale(n_dof);
            
        end

    end
    
    methods (Access = public) % Abstract implementations
        
        function scale_type = getScaleType(this)
            
            scale_type = TrajScale.PROP_SCALE;
            
        end
        
    end
    
    methods (Access = protected) % Abstract implementations

        function sc = calcScaling(this)

            sc = diag( (this.Yg - this.Y0) ./ (this.Ygd - this.Y0d) );
            
        end
        
        % ------------------------------------------
        
        function inv_sc = calcInvScaling(this)
            
            inv_sc = diag( (this.Ygd - this.Y0d) ./ (this.Yg - this.Y0) );
            
        end
        
    end
    
end
