% autorunEpsiConvert.m
% automation to convert .ASCII to .MAT
%
% Nicole Couto adapted from autorunFastCTDConvert.m
% May 2021
% -------------------------------------------------------------------------

instrument = 'fctd';

% Directories for Epsi:
rawDir      = '/Volumes/FCTD_EPSI/RAW';
rawDirAway  = '/Volumes/MOD_data_1/FCTD_EPSI/RAW_0719/raw';
matDir      = '/Volumes/MOD_data_1/FCTD_EPSI/RAW_0719/mat2';
FCTDmatDir  = '/Volumes/MOD_data_1/FCTD_EPSI/RAW_0719/FCTDmat';

% % Directories for Epsi:
% rawDir      = '/Volumes/Berry/blt_epsi/0708_timeseries_testAutoRun/';
% rawDirAway  = '/Volumes/Berry/blt_epsi/0708_timeseries_testAutoRun/raw';
% matDir      = '/Volumes/Berry/blt_epsi/0708_timeseries_testAutoRun/mat';
% FCTDmatDir  = '/Volumes/Berry/blt_epsi/0708_timeseries_testAutoRun/FCTDmat';

% % Directories for FCTD dye chase:
% rawDir      = '/Volumes/FCTD_EPSI/RAW_FCTD';
% rawDirAway  = '/Volumes/MOD_data_1/FCTD_EPSI/RAW_0704/raw';
% matDir      = '/Volumes/MOD_data_1/FCTD_EPSI/RAW_0704/mat';
% FCTDmatDir  = '/Volumes/MOD_data_1/FCTD_EPSI/RAW_0704/FCTDmat';

% % Test rsync only files ending in *.raw
% rawDir = '/Volumes/Public/DY132/SIO-MOD/data/epsi/0624_blt_transect_1/RAW_realtime_0624';
% rawDirAway = '/Volumes/Public/DY132/SIO-MOD/data/epsi/0624_blt_transect_1/RAW_realtime_0624/raw';
% matDir = '/Volumes/Public/DY132/SIO-MOD/data/epsi/0624_blt_transect_1/RAW_realtime_0624/mat';
% FCTDmatDir  = '/Volumes/Public/DY132/SIO-MOD/data/epsi/0624_blt_transect_1/RAW_realtime_0624/fctdmat';

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

% Read configuration data:
% First, try reading configuration data from the
% file. If that doesn't work, try reading from a
% % configuration file. fg                                               
% try
%     setupfile=dir(fullfile(rawDir,'*_raw*'));
%     setup=mod_som_read_setup_from_raw(setupfile(1).name);
% catch
    try
        setupfile=dir(fullfile(rawDir,'*config*'));
        setup=mod_som_read_setup_from_config(setupfile.name);
    catch
        error('mod_som_read_setup failed')
    end
%end

% Create Meta_Data
Meta_Data = fill_meta_data(setup);
Meta_Data = read_MetaProcess(Meta_Data,...
    fullfile(Meta_Data.processpath,'Meta_Data_Process','Meta_Data_Process.txt'));
Meta_Data.rawfileSuffix = '.raw';    

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

EpsiConvert_timer = timer;
EpsiConvert_timer_Stop = false;

EpsiConvert_timer.StartFcn = 'disp(''Conversion of Epsi Data begins now!'');';
switch instrument
    case 'fctd'
EpsiConvert_timer.TimerFcn = [...
    'if EpsiConvert_timer_Stop, '...
    'stop(EpsiConvert_timer); '...
    'delete(EpsiConvert_timer); '...
    'else, '...
    'disp([datestr(now) '': Converting ascii data to mat...!'']); '...
    'try, '...
    'Epsi_MakeMatFromRaw(dirs,Meta_Data,''noGrid'',''doFCTD''); '...
    'catch err, '...
    'disp(err); '...
    'end; '...
    ...%'unix(sprintf(''/usr/bin/rsync -av %s %s'',matDir,matDirAway)); '...
    ...%'unix(sprintf(''/usr/bin/rsync -av %s %s'',pdfDir,pdfDirAway)); '...
    ...%'unix(sprintf(''/usr/bin/rsync -av %s %s'',jpgDir,jpgDirAway)); '...
    ...%'unix(sprintf(''/usr/bin/rsync -av %s %s'',pngDir,pngDirAway)); '...
    'end;'];
case 'epsi'
    EpsiConvert_timer.TimerFcn = [...
    'if EpsiConvert_timer_Stop, '...
    'stop(EpsiConvert_timer); '...
    'delete(EpsiConvert_timer); '...
    'else, '...
    'disp([datestr(now) '': Converting ascii data to mat...!'']); '...
    'try, '...
    'Epsi_MakeMatFromRaw(dirs,Meta_Data,''noGrid'',''doFCTD''); '...
    'catch err, '...
    'disp(err); '...
    'end; '...
    ...%'unix(sprintf(''/usr/bin/rsync -av %s %s'',matDir,matDirAway)); '...
    ...%'unix(sprintf(''/usr/bin/rsync -av %s %s'',pdfDir,pdfDirAway)); '...
    ...%'unix(sprintf(''/usr/bin/rsync -av %s %s'',jpgDir,jpgDirAway)); '...
    ...%'unix(sprintf(''/usr/bin/rsync -av %s %s'',pngDir,pngDirAway)); '...
    'end;'];
    end
EpsiConvert_timer.Period = 10;
EpsiConvert_timer.BusyMode = 'drop';
EpsiConvert_timer.Name = 'EpsiConvert_timer';
EpsiConvert_timer.Tag = 'EpsiConvert_timer';
EpsiConvert_timer.StopFcn = 'clear(''rawDir'',''rawDirAway'',''matDir''); disp([datestr(now) '': Stopped EpsiConvert_timer'']);';
EpsiConvert_timer.ExecutionMode = 'fixedSpacing';
% EpsiConvert_timer.ExecutionMode = 'singleShot';
EpsiConvert_timer.TasksToExecute = Inf;
EpsiConvert_timer.ErrorFcn = 'disp(''%%%%%%%%%%%%% Error %%%%%%%%%%%%%'');';

start(EpsiConvert_timer);