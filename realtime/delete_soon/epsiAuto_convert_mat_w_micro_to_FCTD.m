% To use FCTD_GUI with epsi data, we do things in a slightly different
% order. Instead of automatically making FCTD .mat files when we make .mat
% files as data is streaming in, we just make .mat files. 
%
% Then, in this step, we read the .mat files, calculate microstructure
% variables and THEN make the FCTD .mat files.
%
% Run this timer script in one instance of Matlab, and run FastCTD_GUI in
% another.

%function [] = epsiAuto_convert_raw_to_mat_and_plot(input_struct);

% automation to convert .ASCII to .MAT and plot timeseries in realtime
%
% Nicole Couto adapted from autorunFastCTDConvert.m
% May 2021
% Updated: July 2022
% -------------------------------------------------------------------------
% BEST PRACTICES:
%   Run this as a script inside a parent script that specifies options for
%   the data you wish to process (see EXAMPLE_plot_streaming_data.m)
%
% INPUTS:
% input_struct contains the following fields
%
% Required:
%  .raw_dir = path to directory where the data are streaming in
%  .away_dir = path to directory where data will be copied and where 
%              subdirectories raw, mat, and FCTDmat will be created
%  .Meta_Data_process_file = path to Meta_Data_Process text file
%
% Optional:
%  .str_to_match = String to search for files in raw_dir to copy. Often, we'll 
%                have data from the entire cruise and divide deployments up by 
%                day. You could select to only copy data from Feb 16th by 
%                setting this field to '*EPSI_B_PC2_22_02_16*'; (default '*')
%  .refresh_time_sec = refresh period in seconds (default 5)
%  .version      = version of mod_som_read_epsi_files.m to use (default 4)
%  .starting_dnum = earliest datenum to plot (default 2 weeks before today)
% -------------------------------------------------------------------------
input_struct.data_dir = '/Users/ncouto/Desktop/reprocess_blt_21_0708/deployment/';
input_struct.Meta_Data_process_file = ['/Users/ncouto/GitHub/EPSILOMETER/'...
    'Meta_Data_Process/Meta_Data_Process_blt_2022.txt'];
input_struct.cruise_specifics.blt_2022 = 1;
input_struct.refresh_time_sec = 60*5;

% Check inputs and apply defaults
field_list = fields(input_struct);
% Check that all required fields are defined
if ~any(contains(field_list,'data_dir')) || ~any(contains(field_list,'Meta_Data_process_file'))
  error('You must define ''data_dir,'' and ''Meta_Data_process_file'' as fields in input_struct')
else
  data_dir = input_struct.data_dir;
  Meta_Data_process_file = input_struct.Meta_Data_process_file;
end
% Check for optional inputs and define them
if ~contains(field_list,'str_to_match')
  str_to_match = '*';
else
  str_to_match = input_struct.str_to_match;
end
if ~contains(field_list,'refresh_time_sec')
  refresh_time_sec = 60;
else
  refresh_time_sec = input_struct.refresh_time_sec;
end
if ~contains(field_list,'version ')
  version  = 4;
else
  version  = input_struct.version;
end
if ~contains(field_list,'starting_dnum')
    starting_dnum = now-14;
else
    starting_dnum = input_struct.starting_dnum;
end

% Initialize epsi_class to get Meta_Data
obj = epsi_class(data_dir,input_struct.Meta_Data_process_file);

% Define the directories
dirs.mat       = fullfile(data_dir,'mat');
dirs.fctd_mat  = fullfile(data_dir,'FCTDmat');

% Create FCTD directory if it doesn't exist
if ~exist(dirs.mat,'dir')
    eval([ '!mkdir ' strrep(dirs.mat,' ','\ ')]);
end
if ~exist(dirs.fctd_mat,'dir')
    eval([ '!mkdir ' strrep(dirs.fctd_mat,' ','\ ')]);
end

% Loop through .mat files that have not been converted to FCTDmat. Always
% also include the last FCTDmat file because it may have been made with a
% yet unfinished mat file.
EpsiConvert_timer = timer();
s = false; %Stop switch

EpsiConvert_timer.StartFcn = 'disp(''Begining data conversion!'');';
EpsiConvert_timer.TimerFcn = [...
    'if s, '...
    'stop(EpsiConvert_timer); '...
    'delete(EpsiConvert_timer); '...
    'else, '...
    'disp([datestr(now) '': Adding microstructure to mat files...'']); '...
    'try, '...
    'convert_mat_w_micro_to_FCTD(dirs,obj.Meta_Data,input_struct); '... 
    'catch err, '...
    'display_error_stack(err); '...
    's=1; '...
    'end; '...
    'end;'];
EpsiConvert_timer.Period = refresh_time_sec;
EpsiConvert_timer.BusyMode = 'drop';
EpsiConvert_timer.Name = 'EpsiConvert_timer';
EpsiConvert_timer.Tag = 'EpsiConvert_timer';
EpsiConvert_timer.StopFcn = 'clear(''dirs''); disp([datestr(now) '': Stopped EpsiConvert_timer'']);';
EpsiConvert_timer.ExecutionMode = 'fixedSpacing';
% EpsiConvert_timer.ExecutionMode = 'singleShot';
EpsiConvert_timer.TasksToExecute = Inf;
EpsiConvert_timer.ErrorFcn = 'disp(''%%%%%%%%%%%%% Error %%%%%%%%%%%%%'');';

start(EpsiConvert_timer);


