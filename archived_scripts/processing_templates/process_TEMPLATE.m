% process TEMPLATE mission
%
% Make a copy of this file and edit user choices 1-5 to process your
% epsi deployment
%
% See /Volumes/GoogleDrive/Shared
% drives/MOD-data-Epsilometer/Library/EPSILOMETER/README_Processing_Steps.gdoc
% for more information about Meta_Data and processing steps
% -------------------------------------------------------------------------
%% USER CHOICES (1/5) - Point to the deployment log
% process_dir = directory containing epsi library
% logFile     = path to .csv log file
process_dir='/Volumes/GoogleDrive/Shared drives/MOD-data-Epsilometer/Library/EPSILOMETER/';
logFile = ['/Volumes/GoogleDrive/Shared drives/MOD-data-Epsilometer/epsi/SP2013/Sproul_Drops/d11/Log_d11.csv'];

%% USER CHOICES (2/5) - Choose profiles to process               
% If you only want to process a certain range of profiles, list them here.
% Otherwise, leave commented.
profileList = 1;

%% USER CHOICES (3/5) - Choose processing steps
% 1 = do this step
% 0 = load previously processed data
processingSteps = [0,... %(1)  create Meta_Data
                   0,... %(2)  process raw data
                   0,... %(3)  split data into profiles
                   0,... %(4)  calibrate FP07 for the deployment
                   1,... %(5)  run batch processing
                   0];   %(6)  grid profiles

%% USER CHOICES (4/5) - point to EPSI libraries
addpath(genpath(process_dir))

%% USER CHOICES (5/5) - Define additional processing choices that will be added to Meta_Data

% Visually confirm how profiles are divided?
userConfirmsProfiles = 1;

% Define output paths separate from location of logFile (default is to
% create these directories where logFile is found)
% md.L1path = '/Users/ncouto/Dropbox/SIO/projects/epsi/data/kelp_forest/d4/L1/';
% md.Epsipath = '/Users/ncouto/Dropbox/SIO/projects/epsi/data/kelp_forest/d4/epsi';
% md.CTDpath = '/Users/ncouto/Dropbox/SIO/projects/epsi/data/kelp_forest/d4/ctd';
% md.RAWpath = '/Users/ncouto/Dropbox/SIO/projects/epsi/data/kelp_forest/d4/raw';

md.PROCESS.tscan=3;
md.PROCESS.Fs_ctd=8;
md.PROCESS.ctd_fc=45;  %45 Hz
md.PROCESS.dz=.25;  
md.PROCESS.fc1=5;
md.PROCESS.fc2=35;
md.PROCESS.Prmin_prof = 3;
md.PROCESS.Prcrit_prof = 10;
% Depth range over which to compute coherence between shear and a3 -
% default is 20% to 80% of water column
md.PROCESS.Prmin = @(x) nanmin(x) + 0.2*range(x);
md.PROCESS.Prmax = @(x) nanmin(x) + 0.8*range(x);

% % add auxillary device field
% md.CTD.name = 'SBE49';
% md.CTD.SN   =  '0237';
% md.CTD.cal_file = [process_dir '/SBE49/' md.CTD.SN '.cal'];

md.MAP.adjustTemp = false; %set true to adjust temperature channel electrical noise to the data


%% ------------------------------------------------------------------------
%
%
% *****  USER SHOULD NOT NEED TO CHANGE ANYTHING BELOW HERE  *****
%
%
% -------------------------------------------------------------------------

%% (1) Create Meta_Data
if processingSteps(1)
    
    Meta_Data = create_metadata_from_deployment_log(logFile,process_dir,md);

    % If you're remaking Meta_Data, you'll need to do the temperature
    % calibration again.
    if processingSteps(5)
        processingSteps(4) = 1;
    end
    
else
    if isfield(md,'L1path')
        try
            load(fullfile(md.L1path,'Meta_Data'));
        catch
            fname = dir(fullfile(md.RAWpath,'Meta*.mat'));
            load(fullfile(md.RAWpath,fname));
        end
    else
        path_mission = fileparts(logFile);  
        try
            load(fullfile(path_mission,'L1','Meta_Data'));
        catch
            fname = dir(fullfile(path_mission,'raw','Meta*.mat'));
            load(fullfile(path_mission,'raw',fname));
        end
    end
end


%% (2) Process raw data
if processingSteps(2)
    mod_epsi_read_rawfiles(Meta_Data);
end

%% (3) Divide data into profiles
if processingSteps(3)
    EPSI_create_profiles(Meta_Data,Meta_Data.PROCESS.Prmin_prof,Meta_Data.PROCESS.Prcrit_prof,userConfirmsProfiles)
    load(fullfile(Meta_Data.L1path,['Profiles_', Meta_Data.deployment]))
elseif ~processingSteps(3) && processingSteps(5)
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

% %% (7) define TF to get a qc flag
% if processingSteps(7)
%     % TODO do not forget to change the listfile for loop back to length(listfile)
%     [MS,minepsi,~]=mod_epsilometer_concat_MS(Meta_Data);
%     % compute TF between a1 and
%     [H1,H2,fH]=create_velocity_tranfer_function(MS,minepsi,Meta_Data);
%     % compute qc flags
%     %  fc1 and fc2 define the frequency range integration for qc.
%     mod_epsilometer_add_sh_quality_flag(Meta_Data,H1,H2,fH)
%     [MS,~,~]=mod_epsilometer_concat_MS(Meta_Data);
% end
% 
% %% (8) Plot scans with highest and lowest epsilon and chi values
% 
% %% (9) Plot binned turbulence
% if processingSteps(7)
%     Map   =mod_epsilometer_grid_turbulence(Meta_Data,MS);
%     mod_epsilometer_grid_plot(Map,Meta_Data)
%     mod_epsilometer_binned_epsilon(Meta_Data,MS)
%     % mod_epsilometer_binned_chi(Meta_Data,MS)
% end
