%% GMP namespace
%   Contains general utility functions defined in the "gmp_" namespace
%

classdef gmp_
    
    %% ========================================
    %% ==============   MATH   ================
    %% ========================================    
    methods (Access = public, Static)

        %% Returns the quaternion product between two quaternions.
        %  @param[in] Q1: First quaterion as a 4x1 vector.
        %  @param[in] Q2: Second quaterion as a 4x1 vector.
        %  @return: the quaternion product Q1*Q2.
        %  \note All quaternions must have the form Q = [u; v] where u is
        %  the scalar and v the 3x1 vector part.
        function Q12 = quatProd(Q1, Q2)

            Q1 = Q1(:);
            Q2 = Q2(:);

            n1 = Q1(1);
            e1 = Q1(2:4);

            n2 = Q2(1);
            e2 = Q2(2:4);

            Q12 = zeros(4,1);
            Q12(1) = n1*n2 - e1'*e2;
            Q12(2:4) = n1*e2 + n2*e1 + cross(e1,e2);

            % if Q12(1)<0, Q12=-Q12; end

        end
        
        
        %% Returns the inverse of a quaternion.
        %  @param[in] Q: Quaterion as a 4x1 vector [u; v], with u the scalar and v the vector part.
        %  @return: the inverse quaternion Q^{-1}.
        function invQ = quatInv(Q)
            
            invQ = zeros(size(Q));

            n = (Q'*Q);
            if (n < 1e-16)
                invQ = zeros(4,1);
                return;
            end

            invQ(1) = Q(1);
            invQ(2:4) = -Q(2:4);

            invQ = invQ / n;

        end
        
        
        %% Returns the logarithm map of a unit quaternion.
        %  @param[in] Q: Unit quaterion as a 4x1 vector [u; v], with u the scalar and v the vector part.
        %  @return: the quaternion logarithm log(Q).
        function log_Q = quatLog(Q)

            n = Q(1);
            e = Q(2:4);
            norm_e = norm(e);

            if (norm_e > 1e-16)
                theta = 2*real(atan2(norm_e,n));
                log_Q = theta*e/norm_e;
            else
                log_Q = zeros(size(e));
            end

        end
        
        
        %% Returns the quaternion exponential map.
        %  @param[in] v: 3x1 vector, representing the product k*theta
        %  of a rotation in unit time.
        %  @return: the quaternion exponential exp(log_Q).
        function Q = quatExp(v)

            norm_v_rot = norm(v);
            theta = norm_v_rot;

            Q = zeros(4,1);

            if (norm_v_rot > 1e-16)
            Q(1) = cos(theta/2);
            Q(2:4) = sin(theta/2)*v/norm_v_rot;
            else
              Q = [1 0 0 0]';
            end

        end


        %% Returns the quaternion difference between two quaternions as Q1*(Q2^{-1})
        %  @param[in] Q1: First quaterion as a 4x1 vector.
        %  @param[in] Q2: Second quaterion as a 4x1 vector.
        %  @return: the quaternion difference Q1*(Q2^{-1}).
        %  \note All quaternions must have the form Q = [u; v] where u is
        %  the scalar and v the 3x1 vector part.
        function Qdiff = quatDiff(Q1, Q2)
            
            if ( dot(Q1,Q2) < 0), Q2 = -Q2; end
            Qdiff = gmp_.quatProd(Q1, gmp_.quatInv(Q2));

        end
        
        function v = quatLogDiff(Q1, Q2)

            v = gmp_.quatLog(gmp_.quatDiff(Q1,Q2));

        end
        
        %% Returns derivative of log given the rotational velocity and orientation (expressed w.r.t. the initial orientation)
        %  @param[in] rotVel: Rotational velocity.
        %  @param[in] Q1: Orientation expressed w.r.t. the initial orientation.
        %  @return: Derivative of log.
        function qdot = rotVel_to_qLogDot(rotVel, Q1)
            
            if ( (1-abs(Q1(1))) <= 1e-15)
                % inv_Jn = eye(3,3);
                qdot = rotVel;
            else
                w = Q1(1);
                v = Q1(2:4);
                norm_v = norm(v);
                k = v / norm_v;
                s_th2 = norm_v;
                c_th2 = w;
                th2 = atan2(s_th2, c_th2);
                % k_cross = [0 -k(3) k(2); k(3) 0 -k(1); -k(2) k(1) 0];
                % inv_Jn = k*k' + (eye(3,3) - k*k')*th2*c_th2/s_th2 - th2*k_cross;
                Pk_rotVel = k*dot(k,rotVel); % project rotVel on k
                qdot = Pk_rotVel + (rotVel - Pk_rotVel)*th2*c_th2/s_th2 - th2*cross(k,rotVel);
            end

            % qdot = inv_Jn * rotVel;
            
%             JqQ = gmp_.jacob_qLog_Q(Q1);
%             qdot = 0.5*JqQ * gmp_.quatProd([0; rotVel], Q1);

        end
        
        
        function rotVel = qLogDot_to_rotVel(qdot, Q1)
            
            if ( (1-abs(Q1(1))) <= 1e-15)
                % Jn = eye(3,3);
                rotVel = qdot;
            else
                w = Q1(1);
                v = Q1(2:4);
                norm_v = norm(v);
                k = v / norm_v;
                s_th2 = norm_v;
                c_th2 = w;
                th2 = atan2(s_th2, c_th2);
                % k_cross = [0 -k(3) k(2); k(3) 0 -k(1); -k(2) k(1) 0];
                % Jn = k*k' + (eye(3,3) - k*k')*s_th2*c_th2/th2 + (s_th2^2/th2)*k_cross;
                Pk_qdot = k*dot(k,qdot); % project rotVel on k
                rotVel = Pk_qdot + (qdot - Pk_qdot)*s_th2*c_th2/th2 + (s_th2^2/th2)*cross(k,qdot);
            end
            % rotVel = Jn * qdot;

%             JQq = gmp_.jacob_Q_qLog(Q1);
%             rotVel = 2 * gmp_.quatProd( JQq*qdot, gmp_.quatInv(Q1) );
%             rotVel = rotVel(2:4);

        end
        

        function qddot = rotAccel_to_qLogDDot(rotAccel, rotVel, Q1)
            
            rotVelQ = [0; rotVel];
            rotAccelQ = [0; rotAccel];

            J = gmp_.jacob_qLog_Q(Q1);
            Jdot = gmp_.jacobDot_qLog_Q(Q1, rotVel);

            qddot = 0.5 * (Jdot * gmp_.quatProd(rotVelQ, Q1) + J * gmp_.quatProd( rotAccelQ+0.5*gmp_.quatProd(rotVelQ,rotVelQ), Q1 ) );

        end
        

        function rotAccel = qLogDDot_to_rotAccel(qddot, rotVel, Q1)

            qdot = gmp_.rotVel_to_qLogDot(rotVel, Q1);
            invQ1 = gmp_.quatInv(Q1);
            J = gmp_.jacob_Q_qLog(Q1);
            Jdot = gmp_.jacobDot_Q_qLog(Q1, rotVel);

            rotAccel = 2 * ( gmp_.quatProd( Jdot*qdot+J*qddot, invQ1 ) );
            rotAccel = rotAccel(2:4);

        end
        
        
        %% Returns the Jacobian from the derivative of log to the derivative of Q.
        %  @param[in] Q1: The orientation w.r.t. the initial orientation.
        %  @return: Jacobian.
        function JQq = jacob_Q_qLog(Q1)

            if ( (1-abs(Q1(1))) <= 1e-15)
                JQq = [zeros(1, 3); eye(3,3)];
                return;
            end

            w = Q1(1);
            v = Q1(2:4);
            norm_v = norm(v);
            k = v / norm_v;
            s_th2 = norm_v;
            c_th2 = w;
            th2 = atan2(s_th2, c_th2);
            Pk = k*k';
            I_Pk = eye(3,3) - Pk;

            JQq = zeros(4,3);
            JQq(1,:) = -0.5 * s_th2 * k';
            JQq(2:4,:) = 0.5 * ( I_Pk*s_th2/th2 + c_th2*Pk );

        end
        
        
        %% Returns the Jacobian from the derivative of Q to the derivative of log.
        %  @param[in] Q1: The orientation w.r.t. the initial orientation.
        %  @return: Jacobian.
        function JqQ = jacob_qLog_Q(Q1)
            
            if ( (1-abs(Q1(1))) <= 1e-15)
                JqQ = [zeros(3,1) eye(3,3)];
                return;
            end

            w = Q1(1);
            v = Q1(2:4);
            norm_v = norm(v);
            k = v / norm_v;
            s_th2 = norm_v;
            c_th2 = w;
            th2 = atan2(s_th2, c_th2);
            Pk = k*k';
            I_Pk = eye(3,3) - Pk;

            JqQ = zeros(3,4);
            JqQ(:,1) = -2 * s_th2 * k;
            JqQ(:,2:4) = 2 * ( I_Pk*th2/s_th2 + c_th2*Pk );
%             JqQ(:,1) = 2*k*(th2*c_th2 - s_th2)/s_th2^2;
%             JqQ(:,2:4) = 2*eye(3,3)*th2/s_th2;

        end
       
        
        %% Returns the time derivative of the Jacobian from the derivative of log to the derivative of Q.
        %  @param[in] Q1: The orientation w.r.t. the initial orientation.
        %  @return: Jacobian time derivative.
        function JqQ_dot = jacobDot_qLog_Q(Q1, rotVel)

            qdot = gmp_.rotVel_to_qLogDot(rotVel, Q1);

            if ( (1-abs(Q1(1))) <= 1e-15)
                JqQ_dot = [-qdot/3 zeros(3,3)];
                return;
            end

            w = Q1(1);
            v = Q1(2:4);
            norm_v = norm(v);
            k = v / norm_v;
            s_th2 = norm_v;
            c_th2 = w;
            th2 = atan2(s_th2, c_th2);
            Pk = k*k';
            I_Pk = eye(3,3) - Pk;
            temp = (th2*c_th2-s_th2)/s_th2^2;

%             JqQ_dot = zeros(3,4);
%             JqQ_dot(:,1) = ((-th2/s_th2 - 2*c_th2*temp/s_th2)*Pk + (temp/th2)*(eye(3,3)-Pk))*qdot;
%             JqQ_dot(:,2:4) = -temp*dot(k,qdot)*eye(3,3);
            
            JqQ_dot = zeros(3,4);
            JqQ_dot(:,1) = -(c_th2*Pk + (s_th2/th2)*I_Pk) * qdot;
            JqQ_dot(:,2:4) = dot(k,qdot)*( -temp*I_Pk - s_th2*Pk ) + (c_th2/th2 - 1/s_th2)*( k*(qdot'*I_Pk) + (I_Pk*qdot)*k' );
            
        end
        
        
        %% Returns the time derivative of the Jacobian from the derivative of Q to the derivative of log.
        %  @param[in] Q1: The orientation w.r.t. the initial orientation.
        %  @return: Jacobian time derivative.
        function JQq_dot = jacobDot_Q_qLog(Q1, rotVel)

            qdot = gmp_.rotVel_to_qLogDot(rotVel, Q1);

            if ( (1-abs(Q1(1))) <= 1e-15)
                JQq_dot = [-qdot'/4; zeros(3,3)];
                return;
            end

            w = Q1(1);
            v = Q1(2:4);
            norm_v = norm(v);
            k = v / norm_v;
            s_th2 = norm_v;
            c_th2 = w;
            th2 = atan2(s_th2, c_th2);
            Pk = k*k';
            I_Pk = eye(3,3) - Pk;
            temp = ((th2*c_th2-s_th2)/th2^2);

            JQq_dot = zeros(4,3);
            JQq_dot(1,:) = -0.25 * qdot' * (c_th2*Pk + (s_th2/th2)*I_Pk);
            JQq_dot(2:4,:) = 0.25*dot(k,qdot)*( temp*I_Pk - s_th2*Pk ) + 0.25*temp*( k*(qdot'*I_Pk) + (I_Pk*qdot)*k' );
            
        end
        
    end
    
    %% ======================================
    %% ==============   IO   ================
    %% ======================================
    methods (Access = public, Static)
        
        %% Write the GMP model to a file.
        % @param[in] gmp: Pointer to a @GMP or @GMPo object.
        % @param[in] fid: Filename string or object of type @FileIO associated with the file.
        % @param[in] prefix: The prefix that will be used for writing the names of all GMP params (optional, default="").
        function write(gmp, fid, prefix)
            
            if (nargin < 3), prefix=''; end
            
            c = class(gmp);
            if ( strcmpi(c, 'GMP') ), gmp_.writeGMP(gmp, fid, prefix)
            elseif ( strcmpi(c, 'GMPo') ), gmp_.writeGMPo(gmp, fid, prefix)
            else, error('[gmp_::write]: Unsupported class type');
            end
            
        end
        
        %% Reads the GMP model from a file.
        % @param[in] gmp: Pointer to a @GMP or @GMPo object.
        % @param[in] fid: Filename string or object of type @FileIO associated with the file.
        % @param[in] prefix: The prefix that will be used for reading the names of all GMP params (optional, default="").
        function read(gmp, fid, prefix)
            
            if (nargin < 3), prefix=''; end
            
            c = class(gmp);
            if ( strcmpi(c, 'GMP') ), gmp_.readGMP(gmp, fid, prefix)
            elseif ( strcmpi(c, 'GMPo') ), gmp_.readGMPo(gmp, fid, prefix)
            else, error('[gmp_::read]: Unsupported class type');
            end
            
        end
        
    end
    
    methods (Access = private, Static)
        
        %% Write the GMP model to a file.
        % @param[in] gmp: Pointer to a @GMP object.
        % @param[in] fid: Filename string or object of type @FileIO associated with the file.
        % @param[in] prefix: The prefix that will be used for writing the names of all GMP params (optional, default="").
        function writeGMP(gmp, fid, prefix)
            
            if (nargin < 3), prefix=''; end
            
            if (ischar(fid))
               filename = fid;
               fid = FileIO(filename, bitor(FileIO.out,FileIO.trunc));
            end
            
            N_kernels = gmp.numOfKernels();
            n_dofs = gmp.numOfDoFs();
            
            fid.write([prefix 'weights'], gmp.W);
            fid.write([prefix 'damping'], gmp.D);
            fid.write([prefix 'stiffness'], gmp.K);
            fid.write([prefix 'N_kernels'], uint32(N_kernels) );
            fid.write([prefix 'N_DoFs'], uint32(n_dofs) );
            fid.write([prefix 'scale_type'], int32(gmp.traj_sc.getScaleType()) );
            fid.write([prefix 'c'], gmp.c);
            fid.write([prefix 'h'], gmp.h);
            
        end
        
        %% Reads the GMP model from a file.
        % @param[in] gmp: Pointer to a @GMP object.
        % @param[in] fid: Filename string or object of type @FileIO associated with the file.
        % @param[in] prefix: The prefix that will be used for reading the names of all GMP params (optional, default="").
        function readGMP(gmp, fid, prefix)
            
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
            
            gmp.Y0d = gmp.W*gmp.regressVec(0);
            gmp.Ygd = gmp.W*gmp.regressVec(1);
            
            gmp.setY0(gmp.Y0d);
            gmp.setGoal(gmp.Ygd);
            
            n_dofs = gmp.numOfDoFs();
            gmp.y_dot = zeros(n_dofs,1);
            gmp.z_dot = zeros(n_dofs,1);
            
            if (scale_type == TrajScale.PROP_SCALE), gmp.setScaleMethod(TrajScale_Prop(n_dofs));
            elseif (scale_type == TrajScale.ROT_MIN_SCALE), gmp.setScaleMethod(TrajScale_Rot_min());
            elseif (scale_type == TrajScale.ROT_WB_SCALE), gmp.setScaleMethod(TrajScale_Rot_wb());
            else, error(['[GMP_IO::read]: Unsupported scale type ''' num2str(scale_type) '''...']);
            end
            
        end
        
        %% Write the GMP model to a file.
        % @param[in] gmp: Pointer to a @GMPo object.
        % @param[in] fid: Filename string or object of type @FileIO associated with the file.
        % @param[in] prefix: The prefix that will be used for writing the names of all GMP params (optional, default="").
        function writeGMPo(gmp, fid, prefix)
            
            if (nargin < 3), prefix=''; end
            
            if (ischar(fid))
                filename = fid;
                fid = FileIO(filename, bitor(FileIO.out,FileIO.trunc));
                gmp_.writeGMP(gmp, fid, prefix);
                fid.write([prefix 'Qd0'], gmp.Qd0);
            else
                gmp_.writeGMP(gmp, fid, prefix);
            end
            
        end
        
        %% Reads the GMP model from a file.
        % @param[in] gmp: Pointer to a @GMPo object.
        % @param[in] fid: Filename string or object of type @FileIO associated with the file.
        % @param[in] prefix: The prefix that will be used for reading the names of all GMP params (optional, default="").
        function readGMPo(gmp, fid, prefix)
            
            if (nargin < 3), prefix=''; end
            
            if (ischar(fid))
                filename = fid;
                fid = FileIO(filename, FileIO.in);
                gmp_.readGMP(gmp, fid, prefix);
                gmp.Qd0 = fid.read([prefix 'Qd0']);
            else
                gmp_.readGMP(gmp, fid, prefix);
            end
            
        end
        
    end

end
