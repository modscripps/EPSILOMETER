% process TLC mission
%
% -------------------------------------------------------------------------
%% Define process directory and data directory
% process_dir = directory containing epsi library
% configFile  = path to config file
process_dir = '~/GitHub/EPSILOMETER/';
metaData_processFile = fullfile(process_dir,'Meta_Data_Process','Meta_Data_Process_tlc_2023.txt');
data_dir = '~/Desktop/reprocess_tlc';
deployment_list = {'02_13_day5_epsi_d1';...
    '02_13_day5_epsi_d2';...
    '02_13_day5_epsi_d3'};

% Do you want to grid the data?
gridData = 01;
P = 0:1:400;

% Add processing library
addpath(genpath(process_dir))

for  d=1:length(deployment_list)

    deployment_dir = fullfile(data_dir,deployment_list{d});

    %% Initialize epsi class and read the Meta_Data Process info for BLT
    ec = epsi_class(deployment_dir,metaData_processFile);
    % ec.Meta_Data.PROCESS.tscan = 3;
    % ec.Meta_Data.PROCESS.dz = 0.5;
    % ec.Meta_Data.PROCESS.nfft = round(ec.Meta_Data.PROCESS.nfft/2);
    % ec.Meta_Data.PROCESS.nfftc = round(ec.Meta_Data.PROCESS.nfftc/2);
    % Meta_Data = ec.Meta_Data;
    % save(fullfile(ec.Meta_Data.paths.data,'Meta_Data.mat'),'Meta_Data');

    %% Convert raw to mat
    ec.f_readData;

    %% Compute turbulence variables and grid to P array (optional)
    switch gridData
        case 0
            ec.f_processNewProfiles
        case 1
            ec.f_processNewProfiles('grid',P)
    end

end










