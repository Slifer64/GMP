clc;
close all;
clear;

set_matlab_utils_path();

load('data/recorded_data.mat');

Ts = 0.002;

q0 = Data.JOINT_POS(:,1);
Pos = Data.TOOL_POS;
Quat = Data.TOOL_ORIENT;

n_data = size(Pos,2);

Q_prev = Quat(:,1);
for j=2:n_data
   a =  Quat(:,j)'*Q_prev;
   
   if (a<0), Quat(:,j) = -Quat(:,j); end
   Quat(:,j) = Quat(:,j) / norm(Quat(:,j));
   
   Q_prev = Quat(:,j);
end


Time = (0:(n_data-1)) * Ts;

RotVel = zeros(3, n_data);

for j=1:n_data-1
    RotVel(:,j) = quatLog( quatProd(Quat(:,j+1), quatInv(Quat(:,j)) ) )/Ts;
end

dPos = zeros(3, n_data);
ddPos = zeros(3, n_data);
RotAccel = zeros(3, n_data);
for i=1:3
   dPos(i,:) = [diff(Pos(i,:)) 0]/Ts;
   ddPos(i,:) = [diff(dPos(i,:)) 0]/Ts;
   RotAccel(i,:) = [diff(RotVel(i,:)) 0]/Ts;
end

Data = struct('Time',Time, 'Pos',Pos, 'Vel',dPos, 'Accel',ddPos, 'Quat',Quat, 'RotVel',RotVel, 'RotAccel',RotAccel);

save('kuka_data.mat', 'Data');

fid = fopen('kuka_train_data.bin','w');
io_.write_mat(q0, fid, true);
io_.write_mat(Time, fid, true);
io_.write_mat(Pos, fid, true);
io_.write_mat(dPos, fid, true);
io_.write_mat(ddPos, fid, true);
io_.write_mat(Quat, fid, true);
io_.write_mat(RotVel, fid, true);
io_.write_mat(RotAccel, fid, true);
fclose(fid);

    



