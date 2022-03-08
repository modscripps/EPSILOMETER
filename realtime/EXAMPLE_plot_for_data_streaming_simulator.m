addpath(genpath('~/GitHub/EPSILOMETER'))
rmpath(genpath('~/GitHub/EPSILOMETER/archived_scripts'))

rawDir = '/Users/ncouto/Desktop/reprocess_blt_21_0708/RAW';
awayDir = '/Users/ncouto/Desktop/reprocess_blt_21_0708/deployment/'; %Make sure there's a config file in here if you need one
str_to_match = '*';
raw_file_suffix = '.raw';
Meta_Data_Process_file = ['/Users/ncouto/GitHub/EPSILOMETER/'...
    'Meta_Data_Process/Meta_Data_Process_blt.txt'];
refresh_time_sec = 1;

epsiAuto_convert_raw_to_mat_and_plot