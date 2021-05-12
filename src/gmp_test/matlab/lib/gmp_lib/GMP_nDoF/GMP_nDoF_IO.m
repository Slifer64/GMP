%% N-DoF GMP class
%  Generalized movement primitive.
%

classdef GMP_nDoF_IO < matlab.mixin.Copyable
    
    methods (Static, Access = public)

        %% Write the GMP model to a file.
        % @param[in] gmp: Pointer to a @GMP_nDoF object.
        % @param[in] fid: Filename string or object of type @FileIO associated with the file.
        % @param[in] prefix: The prefix that will be used for writing the names of all GMP params (optional, default="").
        function write(fid, prefix)
            
            if (nargin < 3), prefix=''; end
            
            if (ischar(fid))
               filename = fid;
               fid = FileIO(filename, bitor(FileIO.out,FileIO.trunc));
            end
            
            N_kernels = uint32(gmp.numOfKernels());
            n_dofs = uint32(gmp.numOfDoFs());
            
            fid.write([prefix 'weights'], gmp.W);
            fid.write([prefix 'damping'], gmp.D);
            fid.write([prefix 'stiffness'], gmp.K);
            fid.write([prefix 'N_kernels'], N_kernels);
            fid.write([prefix 'N_DoFs'], n_dofs);
            fid.write([prefix 'scale_type'], gmp.traj_sc.getScaleType());
            fid.write([prefix 'c'], gmp.c);
            fid.write([prefix 'h'], gmp.h);
            
        end
        
        %% Reads the GMP model from a file.
        % @param[in] gmp: Pointer to a @GMP_nDoF object.
        % @param[in] fid: Filename string or object of type @FileIO associated with the file.
        % @param[in] prefix: The prefix that will be used for reading the names of all GMP params (optional, default="").
        function read(fid, prefix)
            
            if (nargin < 3), prefix=''; end
            
            if (ischar(fid))
               filename = fid;
               fid = FileIO(filename, FileIO.in);
            end
            
            gmp.W = fid.read([prefix 'weights']);
            gmp.D = fid.read([prefix 'damping']);
            gmp.K = fid.read([prefix 'stiffness']);
            %N_kernels = fid.read([prefix 'N_kernels']);
            %n_dofs = fid.read([prefix 'N_DoFs']);
            scale_type = fid.read([prefix 'scale_type']);
            gmp.c = fid.read([prefix 'c']);
            gmp.h = fid.read([prefix 'h']);
            
            n_dofs = gmp.numOfDoFs();
            
            gmp.Y0d = gmp.W*gmp.regressVec(0);
            gmp.Ygd = gmp.W*gmp.regressVec(1);
            gmp.Y0 = gmp.Y0d;
            gmp.Yg = gmp.Ygd;
            gmp.y_dot = zeros(n_dofs,1);
            gmp.z_dot = zeros(n_dofs,1);
            
            if (scale_type == TrajScale.PROP_SCALE), gmp.setScaleMethod( TrajScale_Prop(n_dofs) );
            elseif (scale_type == TrajScale.ROT_MIN_SCALE), gmp.setScaleMethod( TrajScale_Rot_min() );
            elseif (scale_type == TrajScale.ROT_WB_SCALE), gmp.setScaleMethod( TrajScale_Rot_wb() );
            else, error(['[GMP_nDoF_IO::read]: Unsupported scale type ''' num2str(scale_type) '''...\n']);
            end

    end
    
end
