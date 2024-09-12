% epsiAuto_timeseries(input_struct);
% This script does the following on a timer:
%   1. epsiProcess_convert_lastN_raw_to_mat = epsiProcess_convert_new_raw_to_mat
%               - converts the last raw file (or the last two if there's
%               not enough data in the last one)
%               - does NOT save the data as a .mat file
%               - the latest mat data are output in the structure 'matData'
%   2. epsiAuto_get_updated_data
%               - finds the newest matData that is not already stored in
%               the strucuture 'obj'. As data are streaming in, the most
%               recent data file is continuously updated. This function grabs
%               the latest data for plotting
%       2.1 If the new data file has less than 30 seconds of data, merge it
%       with the previous one.
%   3. epsiPlot_timeseries
%               - plots the latest timeseries data
%
% EpsiConvert_timer.Period sets the number of seconds for the timer
%
% Nicole Couto adapted from autorunFastCTDConvert.m
% Summer 2021
% Updated May 2023 -  removed saving data. This just converts .raw to .mat
% to plot but does not save any .mat files
% -------------------------------------------------------------------------
% BEST PRACTICES:
%   Run this as a script inside a parent script that specifies options for
%   the data you wish to process (see RUN_epsiAuto_timeseries.m)
%
% INPUTS:
% input_struct contains the following fields
%
% Required:
%  .raw_dir = path to directory where the data are streaming in
%  .Meta_Data_process_file = path to Meta_Data_Process text file
%
% Optional:
%  .refresh_time_sec = refresh period in seconds (default 5)
%  .sec_to_store  = number of seconds of data that will be stored in memory
%  .sec_to_plot   = number of seconds of data that will be plotted in
%                   timeseries (default 30)
%  .version       = version of mod_som_read_epsi_files.m to use (default 4)
%  .starting_dnum = earliest datenum to plot (default 2 weeks before today)
% -------------------------------------------------------------------------

%% Check inputs and apply defaults
field_list = fields(input_struct);

% Check that all required fields are defined
if ~any(contains(field_list,'raw_dir')) || ~any(contains(field_list,'Meta_Data_process_file'))
  error('You must define ''raw_dir'' and ''Meta_Data_process_file'' as fields in the input structure')
else
  raw_dir = input_struct.raw_dir;
  raw_dir_raw = fullfile(input_struct.raw_dir,'raw');
  Meta_Data_process_file = input_struct.Meta_Data_process_file;
end

% Check for optional inputs and define them
if ~contains(field_list,'refresh_time_sec')
  refresh_time_sec = 5;
else
  refresh_time_sec = input_struct.refresh_time_sec;
end
if ~contains(field_list,'sec_to_store')
  sec_to_store = 120;
else
  sec_to_store = input_struct.sec_to_store;
end
if ~contains(field_list,'sec_to_plot')
  sec_to_plot = 30;
else
  sec_to_plot = input_struct.sec_to_plot;
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

%% Initialize epsi_class in away_dir and create blank structures to fill with data
obj = epsi_class(raw_dir,Meta_Data_process_file);
obj = epsiSetup_make_empty_structure(obj,sec_to_store);

% Add Meta_Data from most recent file 
%obj = find_and_load_meta_data(obj,raw_dir); %NC commented out Aug 2024 - I
%think we get meta data in the previous line

field_list = {'epsi','ctd','alt','vnav','gps'};
for iField=1:length(field_list)
    tMax.(field_list{iField}) = starting_dnum;
end
% Create an axes that will be the input for the first call to
% epsiPlot_epsi_ctd_alt_timeseries. All subsequent calls will reuse the set
% of axes created by that function.
ax = axes;
ax = epsiPlot_timeseries(obj,0,ax);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

EpsiAuto_timer = timer;
s = false;
debug_mode = false;

EpsiAuto_timer.StartFcn = 'disp(''Begining data conversion for timeseries plotting!'');';
EpsiAuto_timer.TimerFcn = [...
    'if s, '...
    'stop(EpsiConvert_timer); '...
    'delete(EpsiConvert_timer); '...
    'else, '...
    'disp([datestr(now) '': Converting most recent raw file to mat...'']); '...
    'try, '...
    'matData = epsiProcess_convert_lastN_raw_to_mat(raw_dir_raw,obj.Meta_Data); '...
    'if range(matData.epsi.time_s)<sec_to_plot && exist(''matDataOld''), '...
    'disp([datestr(now) '': Converting TWO most recent raw files to mat...'']); '...
    'matData = epsiProcess_convert_lastN_raw_to_mat(raw_dir_raw,obj.Meta_Data,2); '...
    'end, '...
    'if debug_mode, '...
    'fname = sprintf([''matData_'',datestr(now,''HHMMSS'')]); '...
    'eval([''save '' fname '' matData'']); '...
    'end, '...
    '[obj,tMax] = epsiAuto_get_updated_data(obj,matData,tMax); '...
    'matDataOld = matData; '...
    'if debug_mode, '...
    'fname = sprintf([''obj_tMax_'',datestr(now,''HHMMSS'')]); '...
    'eval([''save '' fname '' obj tMax'']); '...
    'end; '...
    'catch err, '...
    'display_error_stack(err); '...
    's = 1;'...
    'end; '...
    'try, '...
    'ax = epsiPlot_timeseries(obj,0,ax); '...
    'catch err, '...
    'display_error_stack(err); '...
    'tStop = 1;'...
    'end; '...
    'end;'];
EpsiAuto_timer.Period = refresh_time_sec;
EpsiAuto_timer.BusyMode = 'drop';
EpsiAuto_timer.Name = 'EpsiConvert_timer';
EpsiAuto_timer.Tag = 'EpsiConvert_timer';
EpsiAuto_timer.StopFcn = 'clear(''raw_dir''); disp([datestr(now) '': Stopped EpsiConvert_timer'']);';
EpsiAuto_timer.ExecutionMode = 'fixedSpacing';
% EpsiConvert_timer.ExecutionMode = 'singleShot';
EpsiAuto_timer.TasksToExecute = Inf;
EpsiAuto_timer.ErrorFcn = 'disp(''%%%%%%%%%%%%% Error %%%%%%%%%%%%%'');';

start(EpsiAuto_timer);


