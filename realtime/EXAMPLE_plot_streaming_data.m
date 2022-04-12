EPSILOMETER_path = '/Volumes/FCTD Softwares used in BLT 2022/EPSILOMETER';
addpath(genpath(EPSILOMETER_path))
rmpath(genpath(fullfile(EPSILOMETER_path,'archived_scripts')))

rawDir              = '/Volumes/FCTD_EPSI/RAW/'; %Where data is streaming in
awayDir             = '/Volumes/FCTD_EPSI/playground/'; % Where you want to copy data. Make sure there's a config file in here if you need one
str_to_match        = '*'; %Specify filename string to match (ex. 22_04_11* if you only want data from April 11th 2022)
raw_file_suffix     = '.raw';
refresh_time_sec    = 1;
version             = 3;
Meta_Data_Process_file = fullfile(EPSILOMETER_path,...
    'Meta_Data_Process','Meta_Data_Process_blt.txt');

epsiAuto_convert_raw_to_mat_and_plot