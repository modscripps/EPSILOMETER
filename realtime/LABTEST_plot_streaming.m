EPSILOMETER_path = '/Volumes/FCTD Softwares used in BLT 2022/EPSILOMETER';
addpath(genpath(EPSILOMETER_path))
rmpath(genpath(fullfile(EPSILOMETER_path,'archived_scripts')))

rawDir = '/Volumes/FCTD_EPSI/RAW/';
awayDir = '/Volumes/FCTD_EPSI/playground/'; %Make sure there's a config file in here if you need one
str_to_match = '*';
raw_file_suffix = '.raw';
Meta_Data_Process_file = fullfile(EPSILOMETER_path,...
    'Meta_Data_Process','Meta_Data_Process_blt.txt');
refresh_time_sec = 1;

epsiAuto_convert_raw_to_mat_and_plot