addpath(genpath('~/GitHub/EPSILOMETER'))
rmpath(genpath('~/GitHub/EPSILOMETER/archived_scripts'))

% REQUIRED INPUTS
input_struct.raw_dir = '/Users/ncouto/Desktop/reprocess_blt_22_0808/DEV1_RAW';
input_struct.away_dir = '/Users/ncouto/Desktop/reprocess_blt_22_0808/DEV3/Realtime/depl01'; %Make sure there's a config file in here if you need one
input_struct.processing_dir = '/Users/ncouto/Desktop/reprocess_blt_22_0808/DEV3/Processing/depl01';
input_struct.Meta_Data_process_file = ['/Users/ncouto/GitHub/EPSILOMETER/'...
    'Meta_Data_Process/Meta_Data_Process_blt_2022.txt'];

% OPTIONAL INPUTS
input_struct.str_to_match = 'EPSI22_08_08_23';
input_struct.refresh_time_sec = 2;
%input_struct.version = 4;
input_struct.starting_dnum = datenum(2021,7,8);

% FCTD directories
input_struct.FCTDmat = '/Users/ncouto/Desktop/reprocess_blt_22_0808/DEV3/FCTDmat/mat/';
input_struct.FCTDrot = '/Users/ncouto/Desktop/reprocess_blt_22_0808/DEV3/FCTDmat/rot/';

%epsiAuto_convert_raw_to_mat_and_plot(input_struct)
epsiAuto_convert_raw_to_mat_and_plot
