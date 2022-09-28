function Meta_Data=create_metadata_from_deployment_log(filename,process_dir,md)
% Meta_Data=create_metadata_from_deployment_log(filename,process_dir,md)
%
% INPUTS
%   filename    = deployment log file
%   process_dir = epsilometer processing directory
%   md          = option structure of user-defined Meta_Data info
% -------------------------------------------------------------------------

path_mission = fileparts(filename);

% filename='~/ARNAUD/SCRIPPS/CRUISE/2020 07 - Student Cruise RV Sproul/Sproul_drops/d0/Log_d0.csv';
Meta_Data.process=process_dir;
Meta_Data.path_mission=path_mission;

%%

fid=fopen(filename);
header=fgetl(fid);

% read filename
disp(['Read ' filename])
nb_line=1;
while(feof(fid)==0)
    A{nb_line}=fgetl(fid);
    nb_line=nb_line+1;
end
fclose(fid);
disp(['Done reading ' filename])
% done reading filename

% define Meta_Data parameters
mission=strsplit(A{1},',');
Meta_Data.mission=mission{2};

vehicle=strsplit(A{3},',');
Meta_Data.vehicle=vehicle{2};

vehicle_name=strsplit(A{4},',');
Meta_Data.vehicle_name=vehicle_name{2};

deployment=strsplit(A{5},',');
Meta_Data.deployment=deployment{2};

NFFT=strsplit(A{6},',');
Meta_Data.PROCESS.nfft=str2double(NFFT{2});
Meta_Data.PROCESS.nfftc=floor(Meta_Data.PROCESS.nfft/3);

Fs_epsi=strsplit(A{7},',');
Meta_Data.PROCESS.Fs_epsi=str2double(Fs_epsi{2});

s1_SN=strsplit(A{8},',');
Meta_Data.epsi.s1.SN=s1_SN{2};

s2_SN=strsplit(A{9},',');
Meta_Data.epsi.s2.SN=s2_SN{2};


t1_SN=strsplit(A{10},',');
Meta_Data.epsi.t1.SN=t1_SN{2};

t2_SN=strsplit(A{11},',');
Meta_Data.epsi.t2.SN=t2_SN{2};

CTD_name=strsplit(A{12},',');
Meta_Data.aux1.name = CTD_name{2};

CTD_SN=strsplit(A{13},',');
Meta_Data.aux1.SN   = sprintf('%04.0f',str2num(CTD_SN{2})); %Make 4-digits just in case it's not in the log file

channels=strsplit(A{14},',');
channels=channels(2:end-1);
channels{1}=channels{1}(2:end);
channels{end}=channels{end}(1:end-1);
Meta_Data.PROCESS.channels=channels;

% Make 'timeseries' from 'channels'. (Add _g to acceleration channels and
for n=1:numel(Meta_Data.PROCESS.channels)
    if contains(Meta_Data.PROCESS.channels(n),'a')
        Meta_Data.PROCESS.timeseries{n} = sprintf('%s_g',Meta_Data.PROCESS.channels{n});
    elseif contains(Meta_Data.PROCESS.channels(n),{'s','t'})
        Meta_Data.PROCESS.timeseries{n} = sprintf('%s_volt',Meta_Data.PROCESS.channels{n});
    elseif contains(Meta_Data.PROCESS.channels(n),'c')
        Meta_Data.PROCESS.timeseries{n} = sprintf('%s_count',Meta_Data.PROCESS.channels{n});
    end
end

% number of channels
Meta_Data.PROCESS.nb_channels=length(Meta_Data.PROCESS.channels);

% recording mode (STREAMING or SD)
recording_mode=strsplit(A{15},',');
Meta_Data.PROCESS.recording_mode=recording_mode{2};

starttime=strsplit(A{16},',');
Meta_Data.starttime=datenum(starttime{2});


MADRE_rev=strsplit(A{17},',');
Meta_Data.MADRE.rev=MADRE_rev{2};
MADRE_SN=strsplit(A{18},',');
Meta_Data.MADRE.SN=MADRE_SN{2};

MAP_rev=strsplit(A{19},',');
Meta_Data.MAP.rev=MAP_rev{2};

MAP_SN=strsplit(A{20},',');
Meta_Data.MAP.SN=MAP_SN{2};

MAP_temp=strsplit(A{21},',');
Meta_Data.MAP.temperature=MAP_temp{2};
MAP_shear=strsplit(A{22},',');
Meta_Data.MAP.shear=MAP_shear{2};


Firm_version=strsplit(A{23},',');
Meta_Data.Firmware.version=Firm_version{2};
Firm_ADCshear=strsplit(A{24},',');
Meta_Data.Firmware.ADCshear=Firm_ADCshear{2};
Firm_ADCFPO7=strsplit(A{25},',');
Meta_Data.Firmware.ADC_FPO7=Firm_ADCFPO7{2};
Firm_ADCaccell=strsplit(A{26},',');
Meta_Data.Firmware.ADC_accellerometer=Firm_ADCaccell{2};
Firm_ADCcond=strsplit(A{27},',');
Meta_Data.Firmware.ADC_cond=Firm_ADCcond{2};
Firm_ADCfilter=strsplit(A{28},',');
Meta_Data.Firmware.ADCfilter=Firm_ADCfilter{2};
% done define Meta_Data parameters

% define the frequency axes for the spectral computation
[~,Meta_Data.PROCESS.fe] = pwelch(0*(1:Meta_Data.PROCESS.nfft),...
    Meta_Data.PROCESS.nfft,[], ...
    Meta_Data.PROCESS.nfft, ...
    Meta_Data.PROCESS.Fs_epsi,'psd');

%% create mission folders
Meta_Data.L1path   = fullfile(Meta_Data.path_mission,'L1');
Meta_Data.Epsipath = fullfile(Meta_Data.path_mission,'epsi');
Meta_Data.CTDpath  = fullfile(Meta_Data.path_mission,'ctd');
Meta_Data.RAWpath  = fullfile(Meta_Data.path_mission,'raw');
Meta_Data.SDRAWpath  = fullfile(Meta_Data.path_mission,'sd_raw');

%% Add user-defined Meta_Data structure variables if they are specified
if exist('md','var')
    
    % Add extra Meta_Data fields
    if isfield(md,'PROCESS')
        processFields = fields(md.PROCESS);
        for iField=1:length(processFields)
            Meta_Data.PROCESS.(processFields{iField})  =  md.PROCESS.(processFields{iField});
        end
    end
    if isfield(md,'CTD')
        ctdFields = fields(md.CTD);
        for iField=1:length(ctdFields)
            Meta_Data.CTD.(ctdFields{iField}) = md.CTD.(ctdFields{iField});
        end
    end
    if isfield(md,'MAP')
        mapFields = fields(md.MAP);
        for iField=1:length(mapFields)
            Meta_Data.MAP.(mapFields{iField}) = md.MAP.(mapFields{iField});
        end
    end
    
    % Add any other Meta_Data fields you defined in md
    fieldList = fields(md);
    rm = cell2mat(cellfun(@(x) ismember(x, {'PROCESS','CTD','MAP'}), fieldList, 'UniformOutput', 0));
    fieldList = fieldList(~rm);
    for iField=1:length(fieldList)
        Meta_Data.(fieldList{iField}) = md.(fieldList{iField});
    end
    
end

%% Create path, add shear Sv numbers, add CTD calibration,
% save CTD calibration files in ctd and Shear calibration files in epsi
Meta_Data=mod_define_meta_data_log(Meta_Data);

%% Save data
save(fullfile(Meta_Data.RAWpath, ...
    ['Meta_' Meta_Data.mission '_' Meta_Data.deployment '.mat']),'Meta_Data')

save(fullfile(Meta_Data.L1path,'Meta_Data.mat'),'Meta_Data')

