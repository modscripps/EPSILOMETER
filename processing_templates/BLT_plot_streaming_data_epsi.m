% REQUIRED INPUTS
% raw_dir = where data are streaming in
% away_dir = where we copy files to see the timeseries in real time
% processing_dir = where we copy files to process them for microstructure
% profiles, etc.
input_struct.raw_dir = '/Users/Shared/EPSI_PROCESSING/RAW/';
%input_struct.raw_dir = 'mod_admin@192.168.1.158:/Users/Shared/FCTD_EPSI/RAW/';
input_struct.away_dir = ['/Users/Shared/EPSI_PROCESSING/Realtime/'...
                        '0810_sta08_d15_ef1_pc2/']; %CHANGE THIS FOR NEW DEPLOYMENT
input_struct.processing_dir = ['/Users/Shared/EPSI_PROCESSING/Processing/'...
                        '0810_sta08_d15_ef1_pc2/']; %CHANGE THIS FOR NEW DEPLOYMENT
input_struct.Meta_Data_process_file = ...
    ['/Volumes/FCTD Softwares used in BLT 2022/EPSILOMETER_FCTD/'...
    'Meta_Data_Process/Meta_Data_Process_blt_2022.txt'];

% OPTIONAL INPUTS
input_struct.str_to_match = 'EPSI22_08_10_210656'; %CHANGE THIS FOR NEW DEPLOYMENT. (Don't use *, this is going into a strfind)
input_struct.refresh_time_sec = 5;
%input_struct.version = 4;
input_struct.starting_dnum = datenum(2022,7,31);

% FCTD directories
input_struct.FCTDmat = '/Users/Shared/EPSI_PROCESSING/FCTD_mat/mat/';
input_struct.FCTDrot = '/Users/Shared/EPSI_PROCESSING/FCTD_mat/rot/';

epsiAuto_convert_raw_to_mat_and_plot
