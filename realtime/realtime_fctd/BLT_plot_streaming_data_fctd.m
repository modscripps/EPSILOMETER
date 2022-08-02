% REQUIRED INPUTS
input_struct.raw_dir = '/Volumes/FCTD_EPSI/RAW/';
input_struct.away_dir = '/Volumes/FCTD_EPSI/Deployments/Dock/'; %CHANGE THIS FOR NEW DEPLOYMENT
input_struct.Meta_Data_process_file = ...
    ['/Volumes/FCTD Softwares used in BLT 2022/EPSILOMETER_FCTD/'...
    'Meta_Data_Process/Meta_Data_Process_blt_2022.txt'];

% OPTIONAL INPUTS
input_struct.str_to_match = 'EPSI22_07_27*'; %CHANGE THIS FOR NEW DEPLOYMENT
input_struct.refresh_time_sec = 2;
%input_struct.version = 4;
input_struct.starting_dnum = datenum(2022,7,29);

fctdAuto_convert_raw_to_mat_and_plot
