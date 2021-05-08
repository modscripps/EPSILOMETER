% process TEMPLATE mission
%
% Make a copy of this file for each deployment you want to process and
% rename it accordingly (e.g
% process_[mission]_[vehicle_name]_[deployment].m)
%
% I want to stress the importance of creating this Meta_Data structure
% BEFORE deployment. It will help the user to remember details if
% engineering notes are not sufficient. I am talking from experience here.
%
% See /Volumes/GoogleDrive/Shared
% drives/MOD-data-Epsilometer/Library/EPSILOMETER/README_Processing_Steps.gdoc
% for more information about Meta_Data and processing steps
% -------------------------------------------------------------------------
%% Point to the deployment log
% process_dir = directory containing epsi library
% logFile     = path to .csv log file
process_dir='/Volumes/GoogleDrive/Shared drives/MOD-data-Epsilometer/Library/EPSILOMETER/';
logFile = ['/Volumes/GoogleDrive/Shared drives/MOD-data-Epsilometer/epsi/SP2013/Sproul_Drops/d11/Log_d11.csv'];

%% Choose processing steps
% 1 = do this step
% 0 = load previously processed data
processingSteps = [1,... %(1)  create Meta_Data
                   1,... %(2)  process raw data
                   1,... %(3)  split data into profiles
                   1,... %(4)  calibrate FP07 for the deployment
                   1,... %(5)  run batch processing
                   1];   %(6)  grid profiles

% If you only want to process a certain range of profiles, list them here.
% Otherwise, leave commented.

% profileList = 1:20;

%% point to EPSI libraries
addpath(genpath(process_dir))

%% Define additional processing choices that will be added to Meta_Data
%md.L1path = '/path/to/L1output';

md.PROCESS.tscan=6;
md.PROCESS.Fs_ctd=8;
md.PROCESS.ctd_fc=45;  %45 Hz
md.PROCESS.dz=.25;  %45 Hz
md.PROCESS.fc1=5;
md.PROCESS.fc2=35;
md.PROCESS.Prmin_prof = 3;
md.PROCESS.Prcrit_prof = 10;
% this  depth range to compute the coherence.
% this WILL choke when depth is varying.
% TODO change the coherence estimation to use a range in percentage of the
% profile
md.PROCESS.Prmin=20; % 20 m for a 100m profile
md.PROCESS.Prmax=80; % 80 m for a 100 profile

% % add auxillary device field
% md.CTD.name = 'SBE49';
% md.CTD.SN   =  '0237';
% md.CTD.cal_file = [process_dir '/SBE49/' md.CTD.SN '.cal'];

md.MAP.adjustTemp = false; %set true to adjust temperature channel electrical noise to the data

%% (1) Create Meta_Data
if processingSteps(1)
    
    Meta_Data = create_metadata_from_deployment_log(logFile,process_dir);
    
    % Add extra Meta_Data fields
    processFields = fields(md.PROCESS);
    ctdFields = fields(md.CTD);
    mapFields = fields(md.MAP);
    
    for iField=1:length(processFields)
        Meta_Data.PROCESS.(processFields{iField}) = md.PROCESS.(processFields{iField});
    end
    for iField=1:length(ctdFields)
        Meta_Data.CTD.(ctdFields{iField}) = md.CTD.(ctdFields{iField});
    end
    for iField=1:length(mapFields)
        Meta_Data.MAP.(mapFields{iField}) = md.MAP.(mapFields{iField});
    end

    % If you're remaking Meta_Data, you'll need to do the temperature
    % calibration again.
    processingSteps(4) = 1;
    
else
    Meta_Data = load_metadata_from_deployment_log(logFile);
end


%% (2) Process raw data
if processingSteps(2)
    mod_epsi_read_rawfiles(Meta_Data);
end

%% (3) Divide data into profiles
if processingSteps(3)
    EPSI_create_profiles(Meta_Data,Meta_Data.PROCESS.Prmin_prof,Meta_Data.PROCESS.Prcrit_prof)
    load(fullfile(Meta_Data.L1path,['Profiles_', Meta_Data.deployment]))
elseif processingSteps(3) && processingSteps(5)
   load(fullfile(Meta_Data.L1path,['Profiles_', Meta_Data.deployment]))
end

%% (4) Calibrate FP07 for the deployment
if processingSteps(4)

    switch Meta_Data.vehicle
        case 'FISH'
            datachoice = 'datadown';
            idxchoice = 'down';
        case 'WW'
            datachoice = 'dataup';
            idxchoice = 'up';
    end
    
    % For the temperature correction, we want a good chunk of data. Try
    % tscan = 50 or if the profile is too short, find the longest profile
    % in the deployment and set tscan to about 0.5-0.8 times the profile
    % length
    % id profile
    [~,id_profile] = max(cellfun(@(C) length(C), EpsiProfiles.(idxchoice)));
    tscanAlternate = floor(0.8*length(EpsiProfiles.(datachoice){id_profile}.epsitime)/320);
    tscanDefault = 50;
    tscan = min([tscanDefault,tscanAlternate]);


    display=0;
    titleStr = strrep([Meta_Data.mission ' ' Meta_Data.vehicle_name ' ' Meta_Data.deployment],'_','\_');
    
    Meta_Data=mod_epsi_temperature_spectra(Meta_Data, ...
        EpsiProfiles.(datachoice){id_profile}, ...
        CTDProfiles.(datachoice){id_profile},...
        titleStr,id_profile,display,tscan);
end

%% (5) Run batch processing - compute turbulence profiles
if processingSteps(5)
    if exist('profileList','var')
        mod_epsilometer_batch_process(Meta_Data,profileList);
    else
        mod_epsilometer_batch_process(Meta_Data);
    end
end

%% (6) Grid profiles
if processingSteps(6)
    grid_Profiles(Meta_Data);
end

%% (7) define TF to get a qc flag
if processingSteps(7)
    % TODO do not forget to change the listfile for loop back to length(listfile)
    [MS,minepsi,~]=mod_epsilometer_concat_MS(Meta_Data);
    % compute TF between a1 and
    [H1,H2,fH]=create_velocity_tranfer_function(MS,minepsi,Meta_Data);
    % compute qc flags
    %  fc1 and fc2 define the frequency range integration for qc.
    mod_epsilometer_add_sh_quality_flag(Meta_Data,H1,H2,fH)
    [MS,~,~]=mod_epsilometer_concat_MS(Meta_Data);
end

%% (8) Plot scans with highest and lowest epsilon and chi values

%% (9) Plot binned turbulence
if processingSteps(7)
    Map   =mod_epsilometer_grid_turbulence(Meta_Data,MS);
    mod_epsilometer_grid_plot(Map,Meta_Data)
    mod_epsilometer_binned_epsilon(Meta_Data,MS)
    % mod_epsilometer_binned_chi(Meta_Data,MS)
end
