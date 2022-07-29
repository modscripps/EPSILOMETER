% autorunEpsiConvert.m
% automation to convert .ASCII to .MAT
%
% Nicole Couto adapted from autorunFastCTDConvert.m
% May 2021
% -------------------------------------------------------------------------

% rawDir needs to be a directory that ONLY has raw files and bench_config
% in it. It cannot have any subdirectories
rawDir =     '/Volumes/FCTD_EPSI/RAW';
rawDirAway =     '/Volumes/FCTD_EPSI/RAW_BLT2/raw';
matDir =     '/Volumes/FCTD_EPSI/RAW_BLT2/mat';

dirs = {rawDir; rawDirAway; matDir};
    
% Read configuration data:
% First, try reading configuration data from the
% file. If that doesn't work, try reading from a
% configuration file. fg                                               
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
% end

% Create Meta_Data
Meta_Data = fill_meta_data(setup);
Meta_Data = read_MetaProcess(Meta_Data,...
    fullfile(Meta_Data.processpath,'EPSILON','Meta_Data_Process.txt'));
     

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

EpsiConvert_timer = timer;
EpsiConvert_timer_Stop = false;

EpsiConvert_timer.StartFcn = 'disp(''Conversion of Epsi Data begins now!'');';
EpsiConvert_timer.TimerFcn = [...
    'if EpsiConvert_timer_Stop, '...
    'stop(EpsiConvert_timer); '...
    'delete(EpsiConvert_timer); '...
    'else, '...
    'disp([datestr(now) '': Converting ascii data to mat...!'']); '...
    'try, '...
    'Epsi_MakeMatFromRaw(dirs,Meta_Data,''noGrid''); '...
    'catch err, '...
    'disp(err); '...
    'end; '...
    ...%'unix(sprintf(''/usr/bin/rsync -av %s %s'',matDir,matDirAway)); '...
    ...%'unix(sprintf(''/usr/bin/rsync -av %s %s'',pdfDir,pdfDirAway)); '...
    ...%'unix(sprintf(''/usr/bin/rsync -av %s %s'',jpgDir,jpgDirAway)); '...
    ...%'unix(sprintf(''/usr/bin/rsync -av %s %s'',pngDir,pngDirAway)); '...
    'end;'];
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