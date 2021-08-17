% autorunEpsiConvert.m
% automation to convert .ASCII to .MAT
%
% Nicole Couto adapted from autorunFastCTDConvert.m
% May 2021
% -------------------------------------------------------------------------

str_to_match = '*07_24*';
%str_to_match = '*';
% ~~~ Make sure Meta_Data_Process is pointing to the correct text file
% (line 75)


% Directories for Epsi:
rawDir      = '/Volumes/Berry/epsi_raw/Copy_of_EPSI_RAW_0724/raw2';
awayDir     = '/Volumes/Berry/epsi_raw/Copy_of_EPSI_RAW_0724';
rawDirAway  = fullfile(awayDir,'raw');
matDir      = fullfile(awayDir,'mat');
FCTDmatDir  = fullfile(awayDir,'FCTDmat');

% % Directories for Epsi:
% rawDir      = '/Volumes/FCTD_EPSI/RAW';
% awayDir     = '/Volumes/MOD_data_1/FCTD_EPSI/FCTD_RAW_0721';
% rawDirAway  = fullfile(awayDir,'raw');
% matDir      = fullfile(awayDir,'mat');
% FCTDmatDir  = fullfile(awayDir,'FCTDmat');

% % Directories for FCTD dye chase:
% rawDir      = '/Volumes/FCTD_EPSI/RAW_FCTD';
% rawDirAway  = '/Volumes/MOD_data_1/FCTD_EPSI/RAW_0704/raw';
% matDir      = '/Volumes/MOD_data_1/FCTD_EPSI/RAW_0704/mat';
% FCTDmatDir  = '/Volumes/MOD_data_1/FCTD_EPSI/RAW_0704/FCTDmat';

% Create directories if they don't exist
if ~exist(rawDirAway,'dir')
    eval([ '!mkdir ' strrep(rawDirAway,' ','\ ')]);
end
if ~exist(matDir,'dir')
    eval([ '!mkdir ' strrep(matDir,' ','\ ')]);
end
if ~exist(FCTDmatDir,'dir')
    eval([ '!mkdir ' strrep(FCTDmatDir,' ','\ ')]);
end

% Group directories to input for Epsi_MakeMatFromRaw
try
    dirs = {rawDir; rawDirAway; matDir; FCTDmatDir};
catch
    dirs = {rawDir; rawDirAway; matDir};
end

cd(matDir)
cd ..
% Read configuration data:
% First, try reading configuration data from the
% file. If that doesn't work, try reading from a
% % configuration file. fg
% try
%     setupfile=dir(fullfile(rawDir,'*_raw*'));
%     setup=mod_som_read_setup_from_raw(setupfile(1).name);
% catch
try
    setupfile=dir(fullfile(awayDir,'*config*'));
    setup=mod_som_read_setup_from_config(setupfile.name);
catch
    error('mod_som_read_setup failed')
end
%end


% Initialize obj with structures big enough to load at least one Epsi .mat
% file into (epsi, ctd, and alt strucutres)
obj = epsiSetup_make_empty_structure;
obj.plot_properties = epsiSetup_set_plot_properties;
% Create Meta_Data
obj.Meta_Data = epsiSetup_fill_meta_data(setup);
obj.Meta_Data = epsiSetup_read_MetaProcess(obj.Meta_Data,...
    fullfile(obj.Meta_Data.processpath,'Meta_Data_Process','Meta_Data_Process_blt.txt'));
obj.Meta_Data.rawfileSuffix = '.raw';
obj.Meta_Data.MATpath = matDir;

% Choose a starting tMax value for getting new data
tMax.epsi = now-14;
tMax.ctd = now-14;
tMax.alt = now-14;

% Create an axes that will be the input for the first call to
% epsiPlot_epsi_ctd_alt_timeseries. All subsequent calls will reuse the set
% of axes created by that function.
ax = axes;
ax = epsiPlot_epsi_ctd_alt_timeseries(obj,0,ax);
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
    'matData = epsiProcess_convert_new_raw_to_mat(dirs,obj.Meta_Data,''noGrid'',''noSync'',''fileStr'',str_to_match); '...
    '[obj,tMax] = epsiAuto_get_updated_data(obj,matData,tMax); '...
    'catch err, '...
    'disp(err); '...
    'end; '...
    'try, '...
    'ax = epsiPlot_epsi_ctd_alt_timeseries(obj,0,ax); '...
    'catch err, '...
    'disp(err); '...
    'end; '...
    ...%'unix(sprintf(''/usr/bin/rsync -av %s %s'',matDir,matDirAway)); '...
    ...%'unix(sprintf(''/usr/bin/rsync -av %s %s'',pdfDir,pdfDirAway)); '...
    ...%'unix(sprintf(''/usr/bin/rsync -av %s %s'',jpgDir,jpgDirAway)); '...
    ...%'unix(sprintf(''/usr/bin/rsync -av %s %s'',pngDir,pngDirAway)); '...
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