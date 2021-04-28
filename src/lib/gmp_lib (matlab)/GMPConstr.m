%% GMP constraint class
%

classdef GMPConstr < matlab.mixin.Copyable
    
    methods (Access = public)
        
        %% Constructor.
        %  @param[in] t: time instant of constraint
        %  @param[in] value: value of constraint
        %  @param[in] type: type of constraint {'=', '>', '<'}.
        function this = GMPConstr(x, value, type)
                
            if (nargin < 3), return; end
            
            this.x = x;
            this.value = value;
            this.type = type;
            
            if (type~='=' && type~='<' && type~='>')
                error('[GMPConstr::GMPConstr()]: Invalid constraint type: ''%s''',type);
            end
            
        end

    end
    
    methods (Access = public, Static)
     
        function constr = createConstr(x_data, bound_data, bound_type)
            
            n = length(x_data);
            constr = repmat(GMPConstr(0,0,'='), n, 1);
            for j=1:n
                constr(j) = GMPConstr(x_data(j), bound_data(:,j), bound_type);
            end
            
        end
        
    end
    
    properties (Access = public)
     
        x % canonical time instant of constraint
        value % value of constraint
        type % type of constraint {'=', '>', '<'}.
        
    end
    
end
