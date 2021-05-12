clc;
close all;
clear;

addpath('utils/');

Ts = 0.002;

[Time, Y_data, dY_data, ddY_data] = createData(0, 10, Ts, [0; 0; 0], [0.65; 0.7; -0.8], [-0.5; 0.6; 0.5]);

Data = struct('Time',Time, 'Pos',Y_data, 'Vel',dY_data, 'Accel',ddY_data);
save('data/pos_data.mat','Data');

%% =======================================================

q_data = Data.Pos;
Q0 = [1 0 0 0]';

n_data = length(Data.Time);
Q_data = zeros(4,n_data);
vRot_data = zeros(3,n_data);
dvRot_data = zeros(3,n_data);
Q_data(:,1) = math_.quatProd(math_.quatExp(q_data(:,1)), Q0);
for j=2:n_data
    Q_data(:,j) = math_.quatProd(math_.quatExp(q_data(:,j)), Q0); 
    Q1 = Q_data(:,j-1);
    Q2 = Q_data(:,j);
    vRot_data(:,j-1) = math_.quatLog(math_.quatDiff(Q2,Q1)) / Ts;
end
for i=1:3, dvRot_data(i,:) = diff([vRot_data(i,:) 0])/Ts; end

Data = struct('Time',Time, 'Quat',Q_data, 'rotVel',vRot_data, 'rotAccel',dvRot_data);
save('data/orient_data.mat','Data');


%% ===================================================================

function [Timed, yd_data, dyd_data, ddyd_data] = createData(t0, tf, Ts, p0, pf, pm)

    Timed = t0:Ts:tf;
    
    yd_data = get5thOrderPol(p0, pf, Timed);

    tm = [0.6 0.5 0.4]*Timed(end);
    sigma_m = [0.9 1.1 1.2];

    ym = zeros(size(yd_data));
    for i=1:size(ym,1)
       ym(i,:) = pm(i)*exp(-0.5*(Timed-tm(i)).^2/sigma_m(i)^2);
    end
%     ym = pm*exp(-0.5*(Timed-tm).^2/0.9^2);
    
    yd_data = yd_data + ym;
    dyd_data = zeros(size(yd_data));
    ddyd_data = zeros(size(yd_data));
    for i=1:3
        dyd_data(i,:) = [diff(yd_data(i,:)) 0]/Ts; 
        ddyd_data(i,:) = [diff(dyd_data(i,:)) 0]/Ts; 
    end

end
