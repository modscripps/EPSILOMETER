% autorunEpsiMakeProfiles.m
%
% - Start in the data directory containing the config file and raw directory.
% - On BLT: Wait to initialize this until after you've done one full drop to the
% bottom: the temperature calibration will be done using data from that
% drop
% -------------------------------------------------------------------------
% Define process directory and data directory
% process_dir = directory containing epsi library
% configFile  = path to config file
process_dir = '/Volumes/Berry/EPSILOMETER/';
configFile = '/Volumes/Berry/blt_epsi/0703_timeseries/bench_config';
% process_dir = '/Volumes/Public/DY132/SIO-MOD/software/epsi-fctd/EPSILOMETER/';
% configFile = '/Volumes/MOD_data_1-1/FCTD_EPSI/RAW_0721/bench_config';
data_dir = fileparts(configFile);
cd(data_dir)
addpath(genpath(process_dir))

% Standard pressure grid onto which we'll interpolate the data
P = [0:1:2200].';

% Initialize epsi class
ec = epsi_class;

% Pick out profile drops from CTD pressure timeseries
ec.f_getProfileIndices;

% Calibrate temperature
ec.f_calibrateTemperature

%% Initialize the timer to auto-update the realtime plot
tStop                            = false;
MakeProfiles_timer               = timer;
MakeProfiles_timer.Period        = 20*60;
MakeProfiles_timer.BusyMode      = 'queue';
MakeProfiles_timer.Name          = 'autorunEpsiMakeProfiles';
MakeProfiles_timer.Tag           = 'autorunEpsiMakeProfiles';
MakeProfiles_timer.StopFcn       = 'disp([datestr(now) '': Stopped autorunEpsiMakeProfiles'']);';
MakeProfiles_timer.ExecutionMode = 'fixedSpacing'; %'singleShot'
MakeProfiles_timer.TasksToExecute= Inf;
MakeProfiles_timer.ErrorFcn      = 'disp(''%%%%%%%%%%%%% Error %%%%%%%%%%%%%'');';
MakeProfiles_timer.StartFcn      = 'disp(''Making epsi profiles...'');';

MakeProfiles_timer.TimerFcn = [...
    'if tStop, '...
        'stop(MakeProfiles_time); '...
        'delete(MakeProfiles_time); '...
    'else, '...
        'disp([datestr(now)]); '...
        'try, '...
            'ec.f_processNewProfiles(''grid'',P); '...
            'blt_plot_gridded_epsi; '...
        'catch err, '...
             'disp(err), '...
             'for j = 1:length(err.stack), '...
                'disp([num2str(j) '' '' err.stack(j).name '' '' num2str(err.stack(j).line)]); '...
            'end; '...
        'end; '...
    'end;'];

start(MakeProfiles_timer);


%             disp(['So... this is the error for tranlating file ' myASCIIfiles(i).name]);
%             disp(err);
%             for j = 1:length(err.stack)
%                 disp([num2str(j) ' ' err.stack(j).name ' ' num2str(err.stack(j).line)]);
%             end



