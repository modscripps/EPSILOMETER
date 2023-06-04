% RUN_epsiAuto_timeseries.m
%
%  This is a wrapper script to set up input values and run
%  epsiAuto_timeseries.m
%
% CREATE:
% input_struct, a structure containing the following fields:
% 
% Required:
%  .raw_dir = path to directory where the data are streaming in
%  .Meta_Data_process_file = path to Meta_Data_Process text file
%
% Optional:
%  .refresh_time_sec = refresh period in seconds (default 5)
%  .sec_to_store  = number of seconds of data that will be stored in memory
%                  (default 120 seconds)
%  .sec_to_plot   = number of seconds of data that will be plotted in
%                   timeseries (default 60)
%  .version       = version of mod_som_read_epsi_files.m to use (default 4)
%  .starting_dnum = earliest datenum to plot (default 2 weeks before today)
% -------------------------------------------------------------------------

% instrument = 'fctd';
instrument = 'epsi';

% These probably will be the same for the whole cruise
input_struct.raw_dir = '/Users/Shared/EPSI_PROCESSING/Realtime_RAW';
input_struct.Meta_Data_process_file = 'Volumes/Software_TFO2023/EPSILOMETER/Meta_Data_Process/MDP_tfo_2023.txt';
input_struct.refresh_time_sec = 2;

switch instrument
    case 'epsi'
        epsiAuto_timeseries
    case 'fctd'
        fctdAuto_timeseries
end