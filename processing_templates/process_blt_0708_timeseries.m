% process TEMPLATE mission
%
% -------------------------------------------------------------------------
%% Define process directory and data directory
% process_dir = directory containing epsi library
% configFile  = path to config file
process_dir = '/Volumes/Public/DY132/SIO-MOD/software/epsi-fctd/EPSILOMETER/';
configFile = '/Volumes/Berry/blt_epsi/0708_timeseries/bench_config';
data_dir = fileparts(configFile);
cd(data_dir)
addpath(genpath(process_dir))

%% Initialize epsi class
ec = epsi_class;

%% Pick out profile drops from CTD pressure timeseries
% To do: Add this to epsi_class. ec.f_getCtdProfileIndices. 
% Or: Add to subfunction that picks out drops in Epsi_MakeMatFromRaw
makeOrLoad = 'make';

switch makeOrLoad
    case 'make'
        CTD = epsiProcess_make_PressureTimeseries(ec.Meta_Data.MATpath);
    case 'load'
        load(fullfile(ec.Meta_Data.MATpath,'PressureTimeseries.mat'))
        CTD = PressureTimeseries;
        clear PressureTimeseries
end

%% Divide timeseries into profiles
for profNum = 1:length(CTD.startdown)
    try
    profIdx = CTD.startdown(profNum):CTD.enddown(profNum);
    tMin = CTD.ctddnum(profIdx(1));
    tMax = CTD.ctddnum(profIdx(end));
    Profile = ec.f_cropTimeseries(tMin,tMax);
    Profile.profNum = profNum;
    
    saveName = fullfile(data_dir,'profiles',sprintf('Profile%03.0f',profNum));
    eval(['save ' saveName ' Profile'])
    clear Profile
    catch 
    end
end

%% Calibrate FPO7 temperature
% This step finds the longest profile in the deployment and uses it to
% calibrate FPO7 temperature to Seabird temperature. The calibration values
% are saved in Meta_Data.AFE.t1 and .t2
ec.f_calibrateTemperature

% Check calibration values
ec.Meta_Data.AFE.t1
ec.Meta_Data.AFE.t2

%% Compute turbulence variables
for profNum = 1:length(CTD.startdown)
    try
    load(fullfile(ec.Meta_Data.L1path,sprintf('Profile%03.0f',profNum)));
    Profile = ec.f_computeTurbulence(Profile);
    
    saveName = fullfile(data_dir,'profiles',sprintf('Profile%03.0f',profNum));
    eval(['save ' saveName ' Profile'])
    clear Profile   
    catch
    end
end

% Grid profiles
grid_profiles










