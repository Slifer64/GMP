clc;
close all;
clear;

fid = FileIO('data/train_data.bin', FileIO.in);
fid.printHeader();

Data = struct('Time',[], 'Quat',[], 'rotVel',[], 'rotAccel',[]);

Data.Time = fid.read('Timed');
Data.Quat = fid.read('Qd_data');
Data.rotVel = fid.read('vRotd_data');

n_data = length(Data.Time);
Data.rotAccel = zeros(3,n_data);
dt = Data.Time(2) - Data.Time(1);
for i=1:3, Data.rotAccel(i,:) = [diff(Data.rotVel(i,:)) 0]/dt; end

save('data/orient_train_data.mat','Data');

