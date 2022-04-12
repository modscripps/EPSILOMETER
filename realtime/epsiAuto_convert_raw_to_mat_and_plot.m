% autorunEpsiConvert.m
% automation to convert .ASCII to .MAT and plot timeseries in realtime
%
% Nicole Couto adapted from autorunFastCTDConvert.m
% May 2021
% -------------------------------------------------------------------------
% BEST PRACTICES:
%   Run this as a script inside a parent script that specifies options for
%   the data you wish to process (see EXAMPLE_plot_streaming_data.m)
% 
% % Directories for Epsi:
% rawDir      = '/Users/Shared/FCTD_EPSI/RAW'; % rawdir is where the data are
% awayDir     = '/Users/Shared/Beyster_fish2'; % awaydir is copy folder
% 
% str_to_match = '*EPSI_B_PC2_22_02_16_14*';
% Meta_Data_Process_file = fullfile(obj.Meta_Data.paths.process_library,...
%     'Meta_Data_Process','Meta_Data_Process_blt.txt');
% refresh_time_sec = 5;
% -------------------------------------------------------------------------
rawDirAway  = fullfile(awayDir,'raw');
matDir      = fullfile(awayDir,'mat');
FCTDmatDir  = fullfile(awayDir,'FCTDmat');

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
    %dirs = {rawDir; matDir; FCTDmatDir}; %If using 'noSync' in epsiProcess_convert_new_raw_to_mat
catch
    dirs = {rawDir; rawDirAway; matDir};
    %dirs = {rawDir; matDir}; %If using 'noSync' in epsiProcess_convert_new_raw_to_mat
end

cd(matDir)
cd ..
% Read configuration data:
% First, try reading configuration data from the
% file. If that doesn't work, try reading from a
% % configuration file. fg
try %Try reading the first raw file
    setupfile=dir(fullfile(rawDir,'*.raw*'));
    setup=mod_som_read_setup_from_raw(fullfile(setupfile(1).folder,setupfile(1).name));
catch err
    for ii=1:length(err.stack)
        disp(['error in ' err.stack(ii).name ' - line ' num2str(err.stack(ii).line)])
    end
    error('mod_som_read_setup failed')
    
    try %Try reading the first .[suffix] file
        setupfile=dir(fullfile(rawDir,str_to_match,raw_file_suffix));
        setup=mod_som_read_setup_from_raw(setupfile(1).name);
    catch err
        for ii=1:length(err.stack)
            disp(['error in ' err.stack(ii).name ' - line ' num2str(err.stack(ii).line)])
        end
        error('mod_som_read_setup failed')
        
        try %Try reading the config file
            setupfile=dir(fullfile(awayDir,'*config*'));
            setup=mod_som_read_setup_from_config(setupfile.name);
        catch err
            for ii=1:length(err.stack)
                disp(['error in ' err.stack(ii).name ' - line ' num2str(err.stack(ii).line)])
            end
            error('mod_som_read_setup failed')
        end
    end
end


% Initialize obj with structures big enough to load at least one Epsi .mat
% file into (epsi, ctd, and alt strucutres). Here, we're making a strucutre
% that looks like an epsi class so that we can use epsi_class functions.
% It's not an epsi_class, though. It's just a strucutre, so we need to
% manually add some fields that would be automatically added by creating
% an epsi_class.
obj = epsiSetup_make_empty_structure;
obj.plot_properties = epsiSetup_set_plot_properties;
% Create Meta_Data
obj.Meta_Data = epsiSetup_fill_meta_data(setup);
obj.Meta_Data = epsiSetup_read_MetaProcess(obj.Meta_Data,Meta_Data_Process_file);
obj.Meta_Data.rawfileSuffix = raw_file_suffix;
obj.Meta_Data.paths.mat_data = matDir;

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
    'matData = epsiProcess_convert_new_raw_to_mat(dirs,obj.Meta_Data,''doFCTD'',''fileStr'',str_to_match,''version'',version); '...
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
EpsiConvert_timer.Period = refresh_time_sec;
EpsiConvert_timer.BusyMode = 'drop';
EpsiConvert_timer.Name = 'EpsiConvert_timer';
EpsiConvert_timer.Tag = 'EpsiConvert_timer';
EpsiConvert_timer.StopFcn = 'clear(''rawDir'',''rawDirAway'',''matDir''); disp([datestr(now) '': Stopped EpsiConvert_timer'']);';
EpsiConvert_timer.ExecutionMode = 'fixedSpacing';
% EpsiConvert_timer.ExecutionMode = 'singleShot';
EpsiConvert_timer.TasksToExecute = Inf;
EpsiConvert_timer.ErrorFcn = 'disp(''%%%%%%%%%%%%% Error %%%%%%%%%%%%%'');';

start(EpsiConvert_timer);
