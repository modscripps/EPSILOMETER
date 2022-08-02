addpath(genpath('~/GitHub/EPSILOMETER'))
rmpath(genpath('~/GitHub/EPSILOMETER/archived_scripts'))

% REQUIRED INPUTS
input_struct.raw_dir = '/Users/ncouto/Desktop/reprocess_blt_21_0708/RAW';
input_struct.away_dir = '/Users/ncouto/Desktop/reprocess_blt_21_0708/deployment'; %Make sure there's a config file in here if you need one
input_struct.Meta_Data_process_file = ['/Users/ncouto/GitHub/EPSILOMETER/'...
    'Meta_Data_Process/Meta_Data_Process_blt_2021.txt'];

% OPTIONAL INPUTS
input_struct.str_to_match = '_1343';
input_struct.refresh_time_sec = 1;
%input_struct.version = 4;
input_struct.starting_dnum = datenum(2021,7,8);

%epsiAuto_convert_raw_to_mat_and_plot(input_struct)
fctdAuto_convert_raw_to_mat_and_plot
