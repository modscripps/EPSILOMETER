% process TEMPLATE mission
%
% Make a copy of this file and edit user choices 1-5 to process your
% epsi deployment
%
% See /Volumes/GoogleDrive/Sharede
% drives/MOD-data-Epsilometer/Library/EPSILOMETER/README_Processing_Steps.gdoc
% for more information about Meta_Data and processing steps
% -------------------------------------------------------------------------
%% USER CHOICES (1/5) - Point to the deployment log
% process_dir = directory containing epsi library
% configFile  = path to config file
process_dir = '/Volumes/Public/DY132/SIO-MOD/software/epsi-fctd/EPSILOMETER/';
configFile = '/Volumes/GoogleDrive/Shared drives/MOD-data-Epsilometer/epsi/SOMEFE/Cruise/Lab/data/epsi/bench/bltbench/bench_config';

%% USER CHOICES (2/5) - Choose profiles to process               
% If you only want to process a certain range of profiles, list them here.
% Otherwise, leave commented.
%profileList = 1;

%% USER CHOICES (3/5) - Choose processing steps
% 1 = do this step
% 0 = load previously processed data
processingSteps =  [1,... %(1)  create epsiClass including Meta_Data, process raw data
                   1,... %(2)  add user-defined fields to Meta_Data
                   1,... %(3)  add raw epsi and ctd data to epsiClass
                   1,... %(4)  split data into profiles
                   0,... %(5)  calibrate FP07 for the deployment
                   1,... %(6)  run batch processing
                   0];   %(7)  grid profiles

%% USER CHOICES (4/5) - Define additional processing choices that will be added to Meta_Data

% Visually confirm how profiles are divided?
md.PROCESS.userConfirmsProfiles = 1;

% Define output paths separate from location of logFile (default is to
% create these directories where logFile is found)
% md.L1path = '/Users/ncouto/Dropbox/SIO/projects/epsi/data/kelp_forest/d4/L1/';
% md.Epsipath = '/Users/ncouto/Dropbox/SIO/projects/epsi/data/kelp_forest/d4/epsi';
% md.CTDpath = '/Users/ncouto/Dropbox/SIO/projects/epsi/data/kelp_forest/d4/ctd';
% md.RAWpath = '/Users/ncouto/Dropbox/SIO/projects/epsi/data/kelp_forest/d4/raw';

md.PROCESS.tscan        = 6;
md.PROCESS.nfft         = 1024;
md.PROCESS.nfftc        = floor(md.PROCESS.nfft/3);
md.PROCESS.Fs_ctd       = 16;
md.PROCESS.ctd_fc       = 45;  %45 Hz
md.PROCESS.dz           = .25;  
md.PROCESS.fc1          = 5;
md.PROCESS.fc2          = 35;
md.PROCESS.Prmin_prof   = 3;
md.PROCESS.Prcrit_prof  = 10;
% Depth range over which to compute coherence between shear and a3 -
% default is 20% to 80% of water column
md.PROCESS.Prmin = @(x) nanmin(x) + 0.2*range(x);
md.PROCESS.Prmax = @(x) nanmin(x) + 0.8*range(x);

% Temporary, because test data has no probe serial number and no cal
md.epsi.s1.Sv = 45;
md.epsi.s2.Sv = 45;
md.epsi.t1.dTdV = 2.86;
md.epsi.t2.dTdV = 2.86;

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

%% Point to EPSI libraries and change directory to where config file lives
addpath(genpath(process_dir))
path_mission = fileparts(configFile);
cd(path_mission)

%% (1) Create Meta_Data
if processingSteps(1)
    
    ec = epsi_class;

    % If you're remaking Meta_Data, you'll need to do the temperature
    % calibration again.
    if processingSteps(5)
        processingSteps(4) = 1;
    end
    
else
    % For now, all .mat files are saved in same directory. Later, when you
    % add the directory stucutre, uncomment the commented lines
    ec = epsi_class;
    load(fullfile(path_mission,'Meta_Data'));
    ec.Meta_Data = Meta_Data;
    clear Meta_Data
    
%     if isfield(md,'L1path')
%         try
%             load(fullfile(md.L1path,'Meta_Data'));
%         catch
%             fname = dir(fullfile(md.RAWpath,'Meta*.mat'));
%             load(fullfile(md.RAWpath,fname));
%         end
%     else
%         path_mission = fileparts(configFile);  
%         try
%             load(fullfile(path_mission,'L1','Meta_Data'));
%         catch
%             fname = dir(fullfile(path_mission,'raw','Meta*.mat'));
%             load(fullfile(path_mission,'raw',fname));
%         end
%     end
end

%% (2) Add user-defined fields to Meta_Data
if processingSteps(2)
if exist('md','var')
    
    % Add extra Meta_Data fields
    if isfield(md,'PROCESS')
        processFields = fields(md.PROCESS);
        for iField=1:length(processFields)
            ec.Meta_Data.PROCESS.(processFields{iField})  =  md.PROCESS.(processFields{iField});
        end
    end
    if isfield(md,'CTD')
        ctdFields = fields(md.CTD);
        for iField=1:length(ctdFields)
            ec.Meta_Data.CTD.(ctdFields{iField}) = md.CTD.(ctdFields{iField});
        end
    end
    if isfield(md,'MAP')
        mapFields = fields(md.MAP);
        for iField=1:length(mapFields)
            ec.Meta_Data.MAP.(mapFields{iField}) = md.MAP.(mapFields{iField});
        end
    end
    
    % Add any other Meta_Data fields you defined in md
    fieldList = fields(md);
    rm = cell2mat(cellfun(@(x) ismember(x, {'PROCESS','CTD','MAP'}), fieldList, 'UniformOutput', 0));
    fieldList = fieldList(~rm);
    for iField=1:length(fieldList)
        ec.Meta_Data.(fieldList{iField}) = md.(fieldList{iField});
    end
    
end
end

%% (3) Add raw epsi and ctd data to epsi class
if processingSteps(3)
    ec.ctd = ec.f_getCtd;
    ec.epsi = ec.f_getEpsi;
end

%% (4) Divide data into profiles
if processingSteps(4)
    
    %First, check for CTD data - if it doesn't exist, invent it
    % (This is useful if you're processing bench test data)
    load(fullfile(ec.Meta_Data.CTDpath,['ctd_' ec.Meta_Data.deployment '.mat']),'ctd');
    if ~isstruct(ctd) %If there is no CTD data at all
        pressureOnly = 0;
        ec.f_inventCTDdata(pressureOnly)
    elseif range(ctd.P)<10 %If there is CTD data on the bench, not profiling
        pressureOnly = 1;
        ec.f_inventCTDdata(pressureOnly)
    end
    
    ec.f_createProfiles
    load(fullfile(ec.Meta_Data.L1path,['Profiles_', ec.Meta_Data.deployment]))
elseif ~processingSteps(4) && processingSteps(5)
   load(fullfile(ec.Meta_Data.L1path,['Profiles_', ec.Meta_Data.deployment]))
end

%% (5) Calibrate FP07 for the deployment
if processingSteps(5)

    switch ec.Meta_Data.vehicle_name
        case 'FISH'
            datachoice = 'datadown';
            idxchoice = 'down';
        case 'WW'
            datachoice = 'dataup';
            idxchoice = 'up';
        otherwise
            datachoice = 'datadown';
            idxchoice = 'down';
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
    titleStr = strrep([ec.Meta_Data.mission ' ' ec.Meta_Data.vehicle_name ' ' ec.Meta_Data.deployment],'_','\_');
    
    ec.Meta_Data=mod_epsi_temperature_spectra_v2(ec.Meta_Data, ...
        EpsiProfiles.(datachoice){id_profile}, ...
        CTDProfiles.(datachoice){id_profile},...
        titleStr,id_profile,display,tscan);
end

%% (6) Run batch processing - compute turbulence profiles
if processingSteps(6)
    if exist('profileList','var')   
        ec.f_computeTurbulence(profileList);
    else
        ec.f_computeTurbulence;
    end
end

%% (7) Grid profiles
if processingSteps(7)
    grid_Profiles(Meta_Data);
end
