% This script does the following on a timer:
%   1. epsiProcess_convert_new_raw_to_mat
%               - converts all the raw files that don't already have a
%               paired mat file to mat
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
nSec = 2; %Seconds from the end of timeseries to center spectra around
tscan = 4; %Length of scan in seconds

% Directories for Epsi:
% (No sync version. Look at epsiAuto_convert_raw_to_mat_and_plot to see how
% to sync raw files to a new directory first)
rawDir      = '/Users/Shared/FCTD_EPSI/RAW';
matDir      = '/Users/Shared/FCTD_EPSI/RAW/mat';
dirs = {rawDir,matDir};

% Choose a starting tMax value for getting new data
% You need an initial starting point. You will be grabbing and plotting
% all data that came in after this initial value. If you set the clock on
% the SOM prior to starting epsi, the time array will be at datenum so
% choose something like 'now - 1' for TMAX. If you did not set the clock on
% the SOM, the time array will be in seconds since you powered on so choose
% something like 0 for TMAX.
TMAX = 0;
%TMAX = now - 1;

% --- END USER CHOICES ----------------------------------------------------



% Create directories if they don't exist
if ~exist(matDir,'dir')
    eval([ '!mkdir ' strrep(matDir,' ','\ ')]);
end

cd(matDir)
cd ..
% Read configuration data:
% First, try reading configuration data from the
% file. If that doesn't work, try reading from a
% % configuration file. fg
try
    setupfile=dir(fullfile(rawDir,'*.raw*'));
    setup=mod_som_read_setup_from_raw(setupfile(1).name);
catch
try
    setupfile=dir('*config*');
    setup=mod_som_read_setup_from_config(setupfile.name);
catch
    error('mod_som_read_setup failed')
end
end

% Initialize obj with structures big enough to load at least one Epsi .mat
% file into (epsi, ctd, and alt strucutres)
obj = epsiSetup_make_empty_structure;
obj.plot_properties = epsiSetup_set_plot_properties;
% Create Meta_Data
obj.Meta_Data = epsiSetup_fill_meta_data(setup);
obj.Meta_Data = epsiSetup_read_MetaProcess(obj.Meta_Data,...
    fullfile(obj.Meta_Data.paths.process_library,'Meta_Data_Process','Meta_Data_Process_blt.txt'));
obj.Meta_Data.rawfileSuffix = '.raw';
obj.Meta_Data.paths.mat_data = matDir;

% Apply TMAX to structure tMax. Since the instruments sample at different
% rates, these will become slightly different from each other in the loop
% as new data come in.
tMax.epsi = TMAX;
tMax.ctd = TMAX;
tMax.alt = TMAX;

% Create an axes that will be the input for the first call to
% epsiPlot_epsi_ctd_alt_timeseries. All subsequent calls will reuse the set
% of axes created by that function.
ax = axes;
[~,~,~,ax] = epsiPlot_spectra_at_tMid(obj,2,4);
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
            'matData = epsiProcess_convert_new_raw_to_mat(dirs,obj.Meta_Data,''noGrid'',''noSync''); '...
            'if range(matData.epsi.time_s)<30 && exist(''matDataOld''), '...
                'matData = epsiProcess_merge_mat_files(matDataOld,matData); '...
            'end, '...
            '[obj,tMax] = epsiAuto_get_updated_data(obj,matData,tMax); '...
            'matDataOld = matData; '...
        'catch err, '...
            'disp(err); '...
        'end; '...
        'try, '...
            'tMid=nanmax(obj.epsi.time_s)-nanmin(obj.epsi.time_s)-nSec;, '...
            '[~,~,~,ax] = epsiPlot_spectra_at_tMid(obj,tMid,tscan,30,1,0,1,ax);, '...
        'catch err, '...
            'disp(err); '...
        'end; '...
    'end;'];
EpsiConvert_timer.Period = 5;
EpsiConvert_timer.BusyMode = 'drop';
EpsiConvert_timer.Name = 'EpsiConvert_timer';
EpsiConvert_timer.Tag = 'EpsiConvert_timer';
EpsiConvert_timer.StopFcn = 'clear(''rawDir'',''rawDirAway'',''matDir''); disp([datestr(now) '': Stopped EpsiConvert_timer'']);';
EpsiConvert_timer.ExecutionMode = 'fixedSpacing';
% EpsiConvert_timer.ExecutionMode = 'singleShot';
EpsiConvert_timer.TasksToExecute = Inf;
EpsiConvert_timer.ErrorFcn = 'disp(''%%%%%%%%%%%%% Error %%%%%%%%%%%%%'');';

start(EpsiConvert_timer);


