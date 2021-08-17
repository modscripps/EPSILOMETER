% process TEMPLATE mission
%
% -------------------------------------------------------------------------
%% Define process directory and data directory
% process_dir = directory containing epsi library
% configFile  = path to config file
process_dir = '/Volumes/Public/DY132/SIO-MOD/software/epsi-fctd/EPSILOMETER/';
configFile = '/Volumes/Berry/blt_epsi/0710_test/bench_config';
metaData_processFile = '/Volumes/Public/DY132/SIO-MOD/software/epsi-fctd/EPSILOMETER/Meta_Data_Process/Meta_Data_Process_blt.txt';
data_dir = fileparts(configFile);
cd(data_dir)
addpath(genpath(process_dir))

%% Initialize epsi class and read the Meta_Data Process info for BLT
ec = epsi_class;
ec.f_read_MetaProcess(metaData_processFile)

%% Convert raw to mat
ec.f_readData;

%% Pick out profile drops from CTD pressure timeseries
ec.f_getProfileIndices

%% Calibrate FPO7 temperature
% This step finds the longest profile in the deployment and uses it to
% calibrate FPO7 temperature to Seabird temperature. The calibration values
% are saved in Meta_Data.AFE.t1 and .t2
ec.f_calibrateTemperature

% Check calibration values
ec.Meta_Data.AFE.t1
ec.Meta_Data.AFE.t2

%% Compute turbulence variables and grid to P array (optional)
gridData = 0;
P = 1000:1:2200;

switch gridData
    case 0
        ec.f_processNewProfiles
    case 1
        ec.f_processNewProfiles('grid',P)
end










