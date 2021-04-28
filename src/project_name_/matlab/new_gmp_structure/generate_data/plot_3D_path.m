clc;
close all;
clear;

addpath('utils/');

load('data/pos_data.mat','Data');
Pos = Data.Pos;

load('data/orient_data.mat','Data');
Quat = Data.Quat;

figure
plot_3Dpath_with_orientFrames(Pos(:,1:10:end), Quat(:,1:10:end), 'numberOfFrames',6, 'LineWidth',3, 'frameLineWidth',2, 'frameScale',0.3, 'animated',true);
