%% Create folder and meta data structure of an Epsilometer mission
%  The schematic of this structure can be found on confluence
function Meta_Data=mod_define_meta_data_log(Meta_Data)

%% Define main names 
disp('Creating Meta_data from log file and creating the epsi, ctd, L1 and raw folders if not present ')

if ~exist(Meta_Data.L1path,'dir')
    % create path
    eval([ '!mkdir ' strrep(Meta_Data.L1path,' ','\ ')]);
end
if ~exist(Meta_Data.Epsipath,'dir')
    % create path
    eval([ '!mkdir ' strrep(Meta_Data.Epsipath,' ','\ ')]);
end
if ~exist(Meta_Data.CTDpath,'dir')
    % create path
    eval([ '!mkdir ' strrep(Meta_Data.CTDpath,' ','\ ')]);
end
if ~exist(Meta_Data.RAWpath,'dir')
    % create path
    eval([ '!mkdir ' strrep(Meta_Data.RAWpath,' ','\ ')]);
end
if ~exist(Meta_Data.SDRAWpath,'dir')
    % create path
    eval([ '!mkdir ' strrep(Meta_Data.SDRAWpath,' ','\ ')]);
end

%% add path fields
Meta_Data.root     = Meta_Data.path_mission;

%add PROCESS fields
if ~isfield(Meta_Data,'PROCESS')
    warning('missing PROCESS meta_data' )
end

% add MADRE fields
if ~isfield(Meta_Data,'MADRE')
    warning('missing MADRE meta_data' )
end
% add MAP fields
if ~isfield(Meta_Data,'MAP')
    warning('missing MAP meta_data' )
end

% add Firmware fields
if ~isfield(Meta_Data,'Firmware')
    warning('missing MADRE Firmware meta_data' )
end

disp('add ADC config')
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

disp('add shear calibration')
Meta_Data.epsi.shearcal_path=fullfile(Meta_Data.process,'CALIBRATION','SHEAR_PROBES');
Meta_Data.epsi=get_shear_calibration(Meta_Data.epsi);    % Calibration number

disp('add MAP filter')
Meta_Data=get_filters_name_MADRE(Meta_Data);


disp('add CTD calibration')
Meta_Data.CALIpath=fullfile(Meta_Data.process,'CALIBRATION','ELECTRONICS');
Meta_Data.SBECALIpath=fullfile(Meta_Data.process,'SBE49');
original_calfile=fullfile(Meta_Data.SBECALIpath,[Meta_Data.aux1.SN '.cal']);
disp('save CTD calibration in /ctd')
copyfile(original_calfile,Meta_Data.CTDpath,'f')
Meta_Data.SBEcal=get_CalSBE(original_calfile);

disp('save shear calibration in /epsi')
file_calib_shear1=fullfile(Meta_Data.epsi.shearcal_path, ...
                           Meta_Data.epsi.s1.SN, ...
                           ['Calibration_' Meta_Data.epsi.s1.SN '.txt']);
file_calib_shear2=fullfile(Meta_Data.epsi.shearcal_path, ...
                           Meta_Data.epsi.s2.SN, ...
                           ['Calibration_' Meta_Data.epsi.s2.SN '.txt']);

copyfile(original_calfile,Meta_Data.CTDpath,'f')
copyfile(original_calfile,Meta_Data.CTDpath,'f')
copyfile(file_calib_shear1,Meta_Data.Epsipath,'f')
copyfile(file_calib_shear2,Meta_Data.Epsipath,'f')

Meta_Data.name_ctd=['ctd_' Meta_Data.deployment];



%MHA 11/7/2019: comment this out too.
%save(fullfile(Meta_Data.RAWpath, ...
%    ['Meta_' Meta_Data.mission '_' Meta_Data.deployment '.mat']),'Meta_Data')
end