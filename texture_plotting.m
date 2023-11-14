%% Plot pole figures/ODF sections from matflow workflow %%

addpath '/Users/user/Documents/MATLAB/mtex-5.7.0'
% addpath '/Users/user/Documents/MATLAB/mtex-5.8.2'
startup_mtex
clear all
clc
setMTEXpref('FontSize', 18);

%% === load phases: === %

% --- User defined --- %
phase = 'Ti_alpha'; % phase keys in phases.m
mrd_max = 10;
% -------------------- %
[cs] = phases(phase); % returns phase parameters


%% ---  Define path to workflow/result file --- %%
% MUST give exact path to dir containing .hdf5 file
HDF_filepath = '/Users/user/Desktop/iCSF-home/postprocessing/Ti64_alpha_randtext_TEST_2023-11-03-143615/task_4_simulate_volume_element_loading/';
% -------------------- %

% load data from workflow. must specify phase
% quat_data = load_workflow(strcat(HDF_filepath,'workflow.hdf5'), phase);

% load data from damask geom_load.hdf5 file. must specify phase
quat_data = load_result(strcat(HDF_filepath,'geom_load.hdf5'), phase);

%% --- Define path to quaternion.txts --- %%
% uncomment line below for
% path_to_txts = '/Users/user/Desktop/iCSF-home/postprocessing/Ti64_alpha_TEST_2023-11-02-101645/task_4_simulate_volume_element_loading/Ti_alpha_oris/';
% -------------------- %

%% loop over increments...
n_incs = 100;

for inc = 1:1:n_incs
    
    increment = string(inc);
    fprintf("\nInc %s:\n", increment);

%     % Define path to quaternion.txt files...
%     ori_path = strcat(path_to_results, phase, '_inc', increment, '_oris.txt')
%     % Read the quaternions from the txt file
%     fid = fopen(ori_path);
%     quat_data = textscan(fid, '%f%f%f%f', 'HeaderLines', 1, 'CollectOutput', 1);
%     quat_data = quat_data{:};
%     fid = fclose(fid);
    q = quaternion(transpose(quat_data(:, :, inc))); % from HDF5
%     q = quaternion(transpose(data)); % from quat.txt
    
    % Estiamte an ODF from the orientations
    ori = orientation(q, cs);
    % calcdensity method
    odf = calcDensity(ori, 'kernel', deLaValleePoussinKernel, 'halfwidth', 10*degree);
    % elliot et als help
%     psi = calcKernel(ori,'Method','KLCV')
%     odf = calcKernelODF(ori,'kernel',psi);

    % for IPF colouring
    IPFkey = ipfColorKey(odf.CS);
    IPFkey.inversePoleFigureDirection = vector3d.(upper('Y'));
    oriColors = IPFkey.orientation2color(ori);

%     plot pole figures
    figure(1)
    newMtexFigure('figSize', 'normal', 'layout', [1,3])
    plot_PF(ori, odf, phase, oriColors, mrd_max)
    saveas(gcf, strcat(HDF_filepath, 'PF_inc', num2str(inc,'%03.f'), '.png'))

  if contains(phase, "beta")
%   plot odf phi2 slices for beta
    figure(2)
    newMtexFigure('figSize', 'normal')
    odfsection(ori, odf, 45, oriColors, mrd_max)
    saveas(gcf, strcat(HDF_filepath, 'ODF_inc', num2str(inc,'%03.f'), '.png'))
  end
end

fprintf("\n==== ALL DONE ====")