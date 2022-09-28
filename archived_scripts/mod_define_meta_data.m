function Meta_Data=mod_define_meta_data(Meta_Data)
% Meta_Data = mod_define_meta_data(Meta_Data)
%
% Create folder and meta data structure of an Epsilometer mission
%  The schematic of this structure can be found on confluence
%
% Get default values of the Meta_Data strucutre. It requires at least a
% Meta_Data.path_mission
% Meta_Data.mission
% Meta_Data.vehicle_name
% Meta_Data.deployment
% Meta_Data.vehicle
% Meta_Data.process_dir
%
% written by ALB 
% update 04/01/2020


%% create mission folders
mission_folder_L0=fullfile(Meta_Data.path_mission,...
                           Meta_Data.mission);
mission_folder_L1=fullfile(mission_folder_L0,...
                        Meta_Data.vehicle_name);
mission_folder_L2=fullfile(mission_folder_L1,...
                        Meta_Data.deployment);
                    
L1path   = fullfile(mission_folder_L2,'L1');
Epsipath = fullfile(mission_folder_L2,'epsi');
CTDpath  = fullfile(mission_folder_L2,'ctd');
RAWpath  = fullfile(mission_folder_L2,'raw');
SDRAWpath  = fullfile(mission_folder_L2,'sd_raw');


%% add path fields
Meta_Data.root     = mission_folder_L2;
Meta_Data.L1path   = L1path;
Meta_Data.Epsipath = Epsipath;
Meta_Data.CTDpath  = CTDpath;
Meta_Data.RAWpath  = RAWpath;
Meta_Data.SDRAWpath  = SDRAWpath;

%add PROCESS fields
if ~isfield(Meta_Data,'PROCESS')
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
end

% add MADRE fields
if ~isfield(Meta_Data,'Hardware')
    Meta_Data.Hardware.SOM.rev='MADREB.0';
    Meta_Data.Hardware.SOM.SN='0002';
    Meta_Data.Hardware.EFE.rev='MADREB.0';
    Meta_Data.Hardware.EFE.SN='0002';
    Meta_Data.Hardware.EFE.temperature='Tdiff';
    Meta_Data.Hardware.EFE.shear='CAmp1.0';

%% add Firmware fields
if ~isfield(Meta_Data,'Firmware')
    Meta_Data.Firmware.version='MADRE2.1';
    Meta_Data.Firmware.sampling_frequency='320Hz';
    Meta_Data.Firmware.ADCshear='Unipolar';
    Meta_Data.Firmware.ADC_FPO7='Unipolar';
    Meta_Data.Firmware.ADC_accellerometer='Unipolar';
    Meta_Data.Firmware.ADC_cond='count';
    Meta_Data.Firmware.ADCfilter='sinc4';
end


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

Meta_Data.epsi.shearcal_path=fullfile(Meta_Data.process_dir,'CALIBRATION','SHEAR_PROBES');
Meta_Data.epsi=get_shear_calibration(Meta_Data.epsi);    % Calibration number

Meta_Data=get_filters_name_MADRE(Meta_Data);

Meta_Data.CALIpath=fullfile(Meta_Data.process_dir,'CALIBRATION','ELECTRONICS');

Meta_Data.PROCESS.h_freq=mod_epsilometer_get_EFE_filters(Meta_Data,Meta_Data.PROCESS.fe);

end