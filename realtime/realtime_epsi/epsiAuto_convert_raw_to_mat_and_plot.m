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
%   .raw_dir = path to directory where the data are streaming in
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
% Check inputs and apply defaults
field_list = fields(input_struct);
% Check that all required fields are defined
if ~any(contains(field_list,'raw_dir')) || ~any(contains(field_list,'away_dir')) || ~any(contains(field_list,'Meta_Data_process_file'))
  error('You must define ''raw_dir,'' ''away_dir,'' and ''Meta_Data_process_file'' as fields in the input structure')
else
  raw_dir = input_struct.raw_dir;
  away_dir = input_struct.away_dir;
  Meta_Data_process_file = input_struct.Meta_Data_process_file;
end
% Check for optional inputs and define them
if ~contains(field_list,'str_to_match')
  str_to_match = '*';
else
  str_to_match = input_struct.str_to_match;
end
if ~contains(field_list,'refresh_time_sec')
  refresh_time_sec = 5;
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

% Define the directories
dirs.raw_incoming = raw_dir;
dirs.raw_copy  = fullfile(away_dir,'raw');
dirs.mat       = fullfile(away_dir,'mat');
dirs.fctd_mat  = fullfile(away_dir,'FCTDmat');

% Create directories if they don't exist
if ~exist(away_dir,'dir')
    eval([ '!mkdir ' strrep(away_dir,' ','\ ')]);
end
if ~exist(dirs.raw_copy,'dir')
    eval([ '!mkdir ' strrep(dirs.raw_copy,' ','\ ')]);
end
if ~exist(dirs.mat,'dir')
    eval([ '!mkdir ' strrep(dirs.mat,' ','\ ')]);
end
if ~exist(dirs.fctd_mat,'dir')
    eval([ '!mkdir ' strrep(dirs.fctd_mat,' ','\ ')]);
end

% Initialize epsi_class in away_dir and create blank structures to fill
% with data
obj = epsi_class(away_dir,Meta_Data_process_file);
obj = epsiSetup_make_empty_structure(obj);


% % Read configuration data:
% % First, try reading configuration data from the
% % file. If that doesn't work, try reading from a
% % % configuration file.
% Meta_Data.paths.raw_data = dirs.raw_incoming;
% Meta_Data = epsiSetup_read_MetaProcess(Meta_Data,Meta_Data_process_file);
% Meta_Data = epsiSetup_get_raw_suffix(Meta_Data);
% raw_file_suffix = Meta_Data.rawfileSuffix;


% try %Try reading the first .[suffix] file
%     setupfile=dir(fullfile(raw_dir,[str_to_match,raw_file_suffix]));
%     setup=mod_som_read_setup_from_raw(fullfile(setupfile(1).folder,setupfile(1).name));
% catch err
%     display('mod_som_read_setup_from_raw failed, trying mod_som_read_setup_from_config...')
%     try %Try reading the config file
%         setupfile=dir(fullfile(away_dir,'*config*'));
%         setup=mod_som_read_setup_from_config(fullfile(setupfile(1).folder,setupfile(1).name));
%     catch err
%         display_error_stack(err)
%         error('mod_som_read_setup_from_config failed')
%     end
% end



% % Initialize obj with structures big enough to load at least one Epsi .mat
% % file into (epsi, ctd, and alt strucutres). Here, we're making a strucutre
% % that looks like an epsi class so that we can use epsi_class functions.
% % It's not an epsi_class, though. It's just a strucutre, so we need to
% % manually add some fields that would be automatically added by creating
% % an epsi_class.
% obj = epsiSetup_make_empty_structure(obj); %make empty structure (default 5 minutes length)
% obj.plot_properties = epsiSetup_set_plot_properties;
% % Create Meta_Data
% obj.Meta_Data = epsiSetup_fill_meta_data(setup);
% obj.Meta_Data = epsiSetup_read_MetaProcess(obj.Meta_Data,Meta_Data_process_file);
% obj.Meta_Data.rawfileSuffix = raw_file_suffix;
% obj.Meta_Data.paths.mat_data = dirs.mat;

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

EpsiConvert_timer = timer();
s = false; %Stop switch

EpsiConvert_timer.StartFcn = 'disp(''Begining data conversion!'');';
EpsiConvert_timer.TimerFcn = [...
    'if s, '...
    'stop(EpsiConvert_timer); '...
    'delete(EpsiConvert_timer); '...
    'else, '...
    'disp([datestr(now) '': Converting ascii data to mat...'']); '...
    'try, '...
    'matData = epsiProcess_convert_new_raw_to_mat(dirs,obj.Meta_Data,''fileStr'',str_to_match,''version'',version); '...
    '[obj,tMax] = epsiAuto_get_updated_data(obj,matData,tMax); '...
    'catch err, '...
    'display_error_stack(err); '...
    's=1; '...
    'end; '...
    'try, '...
    'ax = epsiPlot_timeseries(obj,0,ax); '...
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
