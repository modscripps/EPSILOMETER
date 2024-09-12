% RUN_epsiAuto_process_data.m
%
% This is a wrapper script to set up input values and run
% epsiAuto_process_data.m
% 
% CREATE:
% input_struct, a structure containing the following fields
%
% Required:
%  .raw_dir = path to directory where the data are streaming in 
%  .process_dir = path to directory where data will be copied and where 
%              subdirectories raw, mat, and FCTDmat will be created
%  .Meta_Data_process_file = path to Meta_Data_Process text file
%
% Optional:
%  .str_to_match = (default '*'), String to search for files in raw_dir to copy. Often, we'll 
%                have data from the entire cruise and divide deployments up by 
%                day. You could select to only copy data from Feb 16th by 
%                setting this field to '*EPSI_B_PC2_22_02_16*'; 
%  .refresh_time_sec = (default 5*60), refresh period in seconds 
%  .version       = (default 4), version of mod_som_read_epsi_files.m to use 
% -------------------------------------------------------------------------

% Change these for each deployment
instrument = 'epsi';
input_struct.process_dir = '~/Desktop/DEV3/processed';
input_struct.str_to_match = 'EPSI22';

% -------------------------------------------------------------------------

% These will probably be the same for the whole cruise
input_struct.raw_dir = '~/Desktop/DEV1/raw'; %Where the .modraw files live
input_struct.Meta_Data_process_file = '/Users/ncouto/GitHub/EPSILOMETER/Meta_Data_Process/MDP_blt_2022.txt';
input_struct.refresh_time_sec =  30;
input_struct.cruise_specifics = 'tfo_2023';
switch instrument
    case 'epsi'
        input_struct.depth_array = 0:2000;
    case 'fctd'
        input_struct.depth_array = 0:2000;
end

% Set command window color
set_window_color('cyan')

% Run the processing script on a timer
switch instrument
    case 'epsi'
        epsiAuto_process_data
    case 'fctd'
        fctdAuto_process_data
end

