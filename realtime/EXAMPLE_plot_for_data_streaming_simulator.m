addpath(genpath('~/GitHub/EPSILOMETER'))
rmpath(genpath('~/GitHub/EPSILOMETER/archived_scripts'))

% REQUIRED INPUTS
input_struct.raw_dir = '/Users/ncouto/Desktop/reprocess_blt_21_0708/RAW';
input_struct.away_dir = '/Users/ncouto/Desktop/reprocess_blt_21_0708/deployment'; %Make sure there's a config file in here if you need one
input_struct.Meta_Data_process_file = ['/Users/ncouto/GitHub/EPSILOMETER/'...
    'Meta_Data_Process/Meta_Data_Process_blt.txt'];

% OPTIONAL INPUTS
input_struct.str_to_match = 'EPSI_B_PC2_21_07_08*';
input_struct.refresh_time_sec = 2;
%input_struct.version = 4;
input_struct.starting_dnum = datenum(2021,7,8);

%epsiAuto_convert_raw_to_mat_and_plot(input_struct)
epsiAuto_convert_raw_to_mat_and_plot
