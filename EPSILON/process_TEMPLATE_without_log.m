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
%% Choose processing steps
% 1 = do this step
% 0 = load previously processed data
processingSteps =    [1,... %(1)  create Meta_Data
                      1,... %(2)  process raw data
                      1,... %(3)  split data into profiles
                      1,... %(4)  calibrate FP07 for the deployment
                      1,... %(5)  run batch processing
                      1,... %(6)  grid profiles
                      1,... %(7)  plot scans with highest/lowest epsilon and chi values
                      1,... %(8)  plot binned turbulence
                      1,... %(9)  run EPSI_grid_turbulence
                      1];   %(10) run checkprofile

% If you only want to process a certain range of profiles, list them here.
% Otherwise, leave commented.
% profileList = 1:20;

%% point to EPSI libraries
process_dir='/Volumes/GoogleDrive/Shared drives/MOD-data-Epsilometer/Library/EPSILOMETER/';
addpath(genpath(process_dir))

%% Meta_Data
Meta_Data.path_mission='/Volumes/GoogleDrive/Shared drives/MOD-data-Epsilometer/epsi/EPSI_TUTO/Cruises/';
Meta_Data.mission='MISSION_FOO';
Meta_Data.vehicle_name='epsifisdh';
Meta_Data.deployment='TEMPLATE';
Meta_Data.vehicle='WW';   % 'WireWalker' or 'FISH'
Meta_Data.process_dir=process_dir;

Meta_Data.PROCESS.nb_channels=8;
Meta_Data.PROCESS.channels={'t1','t2','s1','s2','c','a1','a2','a3'};
Meta_Data.PROCESS.recording_mode='SD';
Meta_Data.PROCESS.tscan=6;
Meta_Data.PROCESS.Fs_epsi=325;
Meta_Data.PROCESS.Fs_ctd=8;
Meta_Data.PROCESS.nfft=Meta_Data.PROCESS.tscan*Meta_Data.PROCESS.Fs_epsi;
Meta_Data.PROCESS.nfftc=floor(Meta_Data.PROCESS.nfft/3);
Meta_Data.PROCESS.ctd_fc=45;  %45 Hz
Meta_Data.PROCESS.dz=.25;  %45 Hz
Meta_Data.PROCESS.fc1=5;
Meta_Data.PROCESS.fc2=35;
Meta_Data.PROCESS.Prmin_prof = 3;
Meta_Data.PROCESS.Prcrit_prof = 10;
% this  depth range to compute the coherence.
% this WILL choke when depth is varying.
% TODO change the coherence estimation to use a range in percentage of the
% profile
Meta_Data.PROCESS.Prmin=20; % 20 m for a 100m profile
Meta_Data.PROCESS.Prmax=80; % 80 m for a 100 profile

[~,Meta_Data.PROCESS.fe] = pwelch(0*(1:Meta_Data.PROCESS.nfft),...
                Meta_Data.PROCESS.nfft,[], ...
                Meta_Data.PROCESS.nfft, ...
                Meta_Data.PROCESS.Fs_epsi,'psd');

% add auxillary device field
Meta_Data.CTD.name = 'SBE49';
Meta_Data.CTD.SN   = '0000';
Meta_Data.CTD.cal_file=[Meta_Data.process_dir '/SBE49/' Meta_Data.CTD.SN '.cal'];

% add channels fields
% Meta_Data.epsi.s1.SN='216'; % serial number;
Meta_Data.epsi.s1.SN='000'; % serial number;
Meta_Data.epsi.s2.SN='000'; % serial number;
Meta_Data.epsi.t1.SN='000'; % serial number;
Meta_Data.epsi.t2.SN='000'; % serial number;

Meta_Data.MAP.rev='MAP.0';
Meta_Data.MAP.SN='0003';
Meta_Data.MAP.temperature='';
Meta_Data.MAP.shear='CAmp1.0';
Meta_Data.MAP.adjustTemp = true; %set true to adjust temperature channel electrical noise to the data


%% (1) Create Meta_Data
if processingSteps(1)
    Meta_Data=mod_define_meta_data(Meta_Data);

    % If you're remaking Meta_Data, you'll need to do the temperature
    % calibration again.
    processingSteps(4) = 1;

else
    Meta_Data = load([Meta_Data.L1path,'Meta_Data']);
end

%% (2) Process raw data
if processingSteps(2)
    mod_epsi_read_rawfiles(Meta_Data);
end

%% (3) Divide data into profiles
if processingSteps(3)
    EPSI_create_profiles(Meta_Data,Meta_Data.Prmin_prof,Meta_Data.Prcrit_prof)
    load(fullfile(Meta_Data.L1path,['Profiles_', Meta_Data.deployment]))
else
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


%% (6) define TF to get a qc flag
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

%% (7) Plot scans with highest and lowest epsilon and chi values

%% (8) Plot binned turbulence
if processingSteps(7)
Map   =mod_epsilometer_grid_turbulence(Meta_Data,MS);
mod_epsilometer_grid_plot(Map,Meta_Data)
mod_epsilometer_binned_epsilon(Meta_Data,MS)
% mod_epsilometer_binned_chi(Meta_Data,MS)
end
