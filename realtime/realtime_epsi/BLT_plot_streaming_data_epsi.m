% REQUIRED INPUTS
input_struct.raw_dir = '/Volumes/FCTD_EPSI/RAW/';
input_struct.away_dir = '/Volumes/FCTD_EPSI/Deployments/0804_d2_ef1_pc2/'; %CHANGE THIS FOR NEW DEPLOYMENT
input_struct.Meta_Data_process_file = ...
    ['/Volumes/FCTD Softwares used in BLT 2022/EPSILOMETER_FCTD/'...
    'Meta_Data_Process/Meta_Data_Process_blt_2022.txt'];

% OPTIONAL INPUTS
input_struct.str_to_match = 'EPSI22_08_04_19'; %CHANGE THIS FOR NEW DEPLOYMENT
input_struct.refresh_time_sec = 1;
%input_struct.version = 4;
input_struct.starting_dnum = datenum(2022,7,27);

% FCTD directories
input_struct.FCTDmat = '/Volumes/FCTD_EPSI/FCTD_MAT/mat/';
input_struct.FCTDrot = '/Volumes/FCTD_EPSI/FCTD_MAT/rot/';

%epsiAuto_convert_raw_to_mat_and_plot(input_struct)
epsiAuto_convert_raw_to_mat_and_plot
