function Meta_Data=create_metadata_from_deployment_log_v2(filename)
% Meta_Data=create_metadata_from_deployment_log(filename,process_dir,md)
%
% INPUTS
%   filename    = deployment log file
%   process_dir = epsilometer processing directory
%   md          = option structure of user-defined Meta_Data info
% -------------------------------------------------------------------------

% If epsilib directory is immediately inside EPSILOMETER directory, the
% next two lines should appropriately define the process directory
% ALB to NC: Pretty smart move :).

spltpath=strsplit(path,':');
epsilib_path=spltpath{~cellfun(@isempty, ...
    cellfun(@(x) ...
    strfind(x,'epsilib'),spltpath, ...
    'UniformOutput',false))};

                           
                           
Meta_Data.paths.process_library=fileparts(epsilib_path);
Meta_Data.paths.data=pwd;
Meta_Data.paths.calibration = fullfile(Meta_Data.paths.process_library,'CALIBRATION','ELECTRONICS');

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

Meta_Data=mod_som_get_shear_probe_calibration_v2(Meta_Data);

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
% Some of the old log files have " " around the channel names
if ~isempty(strfind(channels{1},'"'))
channels{1}=channels{1}(2:end);
end
if ~isempty(strfind(channels{end},'"'))
channels{end}=channels{end}(1:end-1);
end
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

% Add ADC info to each channel structure
Meta_Data.epsi.s1.ADCconf=Meta_Data.Firmware.ADCshear; % serial number;
Meta_Data.epsi.s2.ADCconf=Meta_Data.Firmware.ADCshear; % serial number;
Meta_Data.epsi.t1.ADCconf=Meta_Data.Firmware.ADC_FPO7; % serial number;
Meta_Data.epsi.t2.ADCconf=Meta_Data.Firmware.ADC_FPO7; % serial number;
Meta_Data.epsi.c.ADCconf=Meta_Data.Firmware.ADC_cond; % serial number;
Meta_Data.epsi.a1.ADCconf=Meta_Data.Firmware.ADC_accellerometer; % serial number;
Meta_Data.epsi.a2.ADCconf=Meta_Data.Firmware.ADC_accellerometer; % serial number;
Meta_Data.epsi.a3.ADCconf=Meta_Data.Firmware.ADC_accellerometer; % serial number;

Meta_Data.epsi.s1.ADCfilter=Meta_Data.Firmware.ADCfilter; % serial number;
Meta_Data.epsi.s2.ADCfilter=Meta_Data.Firmware.ADCfilter; % serial number;
Meta_Data.epsi.t1.ADCfilter=Meta_Data.Firmware.ADCfilter; % serial number;
Meta_Data.epsi.t2.ADCfilter=Meta_Data.Firmware.ADCfilter; % serial number;
Meta_Data.epsi.c.ADCfilter=Meta_Data.Firmware.ADCfilter; % serial number;
Meta_Data.epsi.a1.ADCfilter=Meta_Data.Firmware.ADCfilter; % serial number;
Meta_Data.epsi.a2.ADCfilter=Meta_Data.Firmware.ADCfilter; % serial number;
Meta_Data.epsi.a3.ADCfilter=Meta_Data.Firmware.ADCfilter; % serial number;


% define the frequency axes for the spectral computation
[~,Meta_Data.PROCESS.fe] = pwelch(0*(1:Meta_Data.PROCESS.nfft),...
    Meta_Data.PROCESS.nfft,[], ...
    Meta_Data.PROCESS.nfft, ...
    Meta_Data.PROCESS.Fs_epsi,'psd');

buildname= @(x,y) fullfile(x,[y '.cal']);
Meta_Data.aux1.CALpath   = fullfile(Meta_Data.paths.process_library,'CALIBRATION','SBE49');
Meta_Data.aux1.CALfile=buildname(Meta_Data.aux1.CALpath,Meta_Data.aux1.SN);
Meta_Data.aux1.cal=get_CalSBE(Meta_Data.aux1.CALfile);

fprintf('Saving Meta_Data in datapath \n')
save(fullfile(Meta_Data.paths.data,'Meta_Data.mat'),'Meta_Data');