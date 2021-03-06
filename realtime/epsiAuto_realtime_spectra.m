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
%   3. epsiPlot_spectra_at_tMid
%               - plots the latest 30 seconds of epsi channel output
%               (t1,t2,s1,s2,a1,a2,a3), the latest 30 seconds of dPdt and
%               frequency spectra of the epsi channels centered on 'nSec'
%               seconds from the end 
% EpsiConvert_timer.Period sets the number of seconds for the timer
%
% Nicole Couto adapted from autorunFastCTDConvert.m
% Summer 2021

% --- USER CHOICES --------------------------------------------------------
totalSec = 120; %Number of seconds of data that will be stored in memory
plotSec = 30; %Number of seconds of data that will be plotted in timeseries
tscan = 4; %Length of scan in seconds
centerScan = tscan/2; %Plotted spectra will be centered tscan/2 seconds from the end of the timeseries
str_to_match = '*';

% Directory containing streaming raw data
rawDir = '/Users/Shared/SOM_APEX/RAW/';
dirs.raw_incoming = '/Users/Shared/SOM_APEX/RAW/';

% Choose time units
time_units = 'dnum'; %uncomment this if you set the datetime on SOM
%time_units = 'seconds'; %uncomment this if you did not set the datetime on SOM

% --- END USER CHOICES ----------------------------------------------------

switch time_units
    case 'seconds'
        TMAX = 0;
    case 'dnum'
        TMAX = datenum(2021,1,1); %Will plot data after this starting point
end

% Read configuration data:
% First, try reading configuration data from the
% file. If that doesn't work, try reading from a
% % configuration file.
Meta_Data.paths.raw_data = dirs.raw_incoming;
[Meta_Data] = epsiSetup_get_raw_suffix(Meta_Data);
raw_file_suffix = Meta_Data.rawfileSuffix;

try %Try reading the first .[suffix] file
    setupfile=dir(fullfile(rawDir,[str_to_match,raw_file_suffix]));
    setup=mod_som_read_setup_from_raw(fullfile(setupfile(1).folder,setupfile(1).name));
catch err
    display_error_stack(err)
    error('mod_som_read_setup failed')
    try %Try reading the config file
        setupfile=dir(fullfile(awayDir,'*config*'));
        setup=mod_som_read_setup_from_config(fullfile(setupfile(1).folder,setupfile(1).name));
    catch err
        display_error_stack(err)
        error('mod_som_read_setup failed')
    end
end

%end

% Initialize obj with structures big enough to load at least one Epsi .mat
% file into (epsi, ctd, and alt strucutres)
obj = epsiSetup_make_empty_structure(totalSec);
obj.plot_properties = epsiSetup_set_plot_properties;
% Create Meta_Data
obj.Meta_Data = epsiSetup_fill_meta_data(setup);
obj.Meta_Data = epsiSetup_read_MetaProcess(obj.Meta_Data,...
    fullfile(obj.Meta_Data.paths.process_library,'Meta_Data_Process','Meta_Data_Process_blt.txt'));
obj.Meta_Data.rawfileSuffix = raw_file_suffix;

% Apply TMAX to structure tMax. Since the instruments sample at different
% rates, these will become slightly different from each other in the loop
% as new data come in.
field_list = {'epsi','ctd','alt','vnav','gps'};
for iField=1:length(field_list)
tMax.(field_list{iField}) = TMAX;
end

% Create an axes that will be the input for the first call to
% epsiPlot_epsi_ctd_alt_timeseries. All subsequent calls will reuse the set
% of axes created by that function.
ax = axes;
[~,~,~,ax] = epsiPlot_spectra_at_tMid(obj,centerScan,tscan);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

EpsiConvert_timer = timer;
tStop = false;

EpsiConvert_timer.StartFcn = 'disp(''Conversion of Epsi Data begins now!'');';
EpsiConvert_timer.TimerFcn = [...
    'if tStop, '...
        'stop(EpsiConvert_timer); '...
        'delete(EpsiConvert_timer); '...
    'else, '...
        'disp([datestr(now) '': Converting ascii data to mat...!'']); '...
        'try, '...
            'matData = epsiProcess_convert_lastN_raw_to_mat(rawDir,obj.Meta_Data); '...
            'if range(matData.epsi.time_s)<plotSec && exist(''matDataOld''), '...
                'matData = epsiProcess_convert_lastN_raw_to_mat(rawDir,obj.Meta_Data,2); '...
            'end, '...
            '[obj,tMax] = epsiAuto_get_updated_data(obj,matData,tMax); '...
            'matDataOld = matData; '...
        'catch err, '...
            'display_error_stack(err); '...
            'tStop = 1;'...
        'end; '...
        'try, '...
            'tMid=nanmax(obj.epsi.time_s)-centerScan;, '...
            '[~,~,~,ax] = epsiPlot_spectra_at_tMid(obj,tMid,tscan,plotSec,1,0,1,ax);, '...
        'catch err, '...
            'display_error_stack(err); '...
            'tStop = 1;'...
        'end; '...
    'end;'];
EpsiConvert_timer.Period = 1;
EpsiConvert_timer.BusyMode = 'drop';
EpsiConvert_timer.Name = 'EpsiConvert_timer';
EpsiConvert_timer.Tag = 'EpsiConvert_timer';
EpsiConvert_timer.StopFcn = 'disp([datestr(now) '': Stopped EpsiConvert_timer'']);';
EpsiConvert_timer.ExecutionMode = 'fixedSpacing';
% EpsiConvert_timer.ExecutionMode = 'singleShot';
EpsiConvert_timer.TasksToExecute = Inf;
EpsiConvert_timer.ErrorFcn = 'disp(''%%%%%%%%%%%%% Error %%%%%%%%%%%%%'');';

start(EpsiConvert_timer);


