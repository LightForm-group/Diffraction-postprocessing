addpath '/Users/user/Documents/MATLAB/mtex-5.7.0'
startup_mtex
clear all
clc

setMTEXpref('FontSize', 32);
%% === load phases: === %

% --- User defined --- %
phase = 'Ti-beta';                                                       % change to 'Ti64_alpha' or 'Ti64_beta'
mrd_max = 2;
% -------------------- %
[cs] = phases(phase);                   % function returns phase parameters
%% Define path to quaternion.txt files...

% --- User defined --- % 
path_to_odf = '/Users/user/Desktop/Diamond_2022/034_Ti64_TIFUN-T4_TD_Deform_910C_1mms-1/fourier-peak-analysis-texture/beta/034_Ti64_TIFUN-T4_TD_Deform_910C_1mms-1_beta_ODF.mat'; % define path
% -------------------- %
%% plot PF...

% import odf
load(path_to_odf);

%     plot pole figure
newMtexFigure('figSize', 'normal', 'layout', [1,4])
plot_PF(odf, phase, cs, mrd_max)

%   plot odf phi2 slices for beta
%     newMtexFigure('figSize', 'huge')
%     odfsection(odf, 45, mrd_max)


% saveas(gcf, texture.png')

disp("==== ALL DONE ====")