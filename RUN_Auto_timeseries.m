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

% -------------------------------------------------------------------------
% --- USER CHOICES --------------------------------------------------------
instrument = 'fctd';
% instrument = 'fctd_tridente';
%instrument = 'epsi';

% Also plot spectra?
include_spectra = 0;

% Meta_Data process file (make sure this file has the correct serial
% numbers for CTD, s1, s2, t1, t2. If you're running fctd, you can leave
% s1, s2, t1, t2 = '115')
Meta_Data_process_file = 'MDP_tfo2024_fctd_0807.txt';

% These probably will be the same for the whole cruise

input_struct.raw_dir = '/Users/ncouto/Desktop/OCEAN';
Meta_Data_process_dir = '~/GitHub/EPSILOMETER/Meta_Data_Process';

%input_struct.raw_dir = '/Users/Shared/EPSI_PROCESSING/TFO2024/Realtime_RAW/';
%Meta_Data_process_dir = '/Volumes/Software_TFO2024/EPSILOMETER/Meta_Data_Process/';
input_struct.Meta_Data_process_file = fullfile(Meta_Data_process_dir,Meta_Data_process_file);
input_struct.refresh_time_sec = 2;

% --- END USER CHOICES ----------------------------------------------------
% -------------------------------------------------------------------------

% Set command window color
set_window_color('yellow')

% Run the realtime plotting script on a timer
switch instrument
    case 'epsi'
        if ~include_spectra
            epsiAuto_timeseries
        elseif include_spectra
            epsiAuto_timeseries_spectra
        end
    case 'fctd'
        fctdAuto_timeseries
    case 'fctd_tridente'
        fctdAuto_timeseries_tridente
end