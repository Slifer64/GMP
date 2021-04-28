%% N-DoF GMP class
%  Generalized movement primitive.
%

classdef GMP_nDoF_IO < matlab.mixin.Copyable
    
    methods (Access = public)
        
        %% GMP constructor.
        %  @param[in] gmp: n_DoF dmp.
        function this = GMP_nDoF_IO(gmp)
                
            this.gmp = gmp;
            
        end

        %% Write the GMP model to a file.
        % @param[in] fid: Filename string or object of type @FileIO associated with the file.
        % @param[in] prefix: The prefix that will be used for writing the names of all GMP params (optional, default="").
        function write(this, fid, prefix)
            
            if (nargin < 3), prefix=''; end
            
            if (ischar(fid))
               filename = fid;
               fid = FileIO(filename, bitor(FileIO.out,FileIO.trunc));
            end
            
            N_kernels = this.gmp.numOfKernels();
            n_dofs = this.gmp.numOfDoFs();
            
            fid.write([prefix 'weights'], this.gmp.W);
            fid.write([prefix 'damping'], this.gmp.D);
            fid.write([prefix 'stiffness'], this.gmp.K);
            fid.write([prefix 'N_kernels'], N_kernels);
            fid.write([prefix 'N_DoFs'], n_dofs);
            fid.write([prefix 'c'], this.gmp.c);
            fid.write([prefix 'h'], this.gmp.h);
            
        end
        
        %% Reads the GMP model from a file.
        % @param[in] fid: Filename string or object of type @FileIO associated with the file.
        % @param[in] prefix: The prefix that will be used for reading the names of all GMP params (optional, default="").
        function read(this, fid, prefix)
            
            if (nargin < 3), prefix=''; end
            
            if (ischar(fid))
               filename = fid;
               fid = FileIO(filename, FileIO.in);
            end
            
            this.gmp.W = fid.read([prefix 'weights']);
            this.gmp.D = fid.read([prefix 'damping']);
            this.gmp.K = fid.read([prefix 'stiffness']);
            %N_kernels = fid.read([prefix 'N_kernels']);
            %n_dofs = fid.read([prefix 'N_DoFs']);
            this.gmp.c = fid.read([prefix 'c']);
            this.gmp.h = fid.read([prefix 'h']);
            
            this.gmp.Y0d = this.gmp.W*this.gmp.regressVec(0);
            this.gmp.Ygd = this.gmp.W*this.gmp.regressVec(1);
            
            this.gmp.setY0(this.gmp.Y0d);
            this.gmp.setGoal(this.gmp.Ygd);
            
            n_dofs = this.gmp.numOfDoFs();
            this.gmp.y_dot = zeros(n_dofs,1);
            this.gmp.z_dot = zeros(n_dofs,1);
            
        end

    end
    
    properties (Access = private)
        
        gmp
        
    end
    
end
