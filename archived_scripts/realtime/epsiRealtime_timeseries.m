% Run this script from Meta_Data.MATpath

instrument = 'epsi';

% Initialize obj with structures big enough to load at least one Epsi .mat
% file into (epsi, ctd, and alt strucutres)
obj = make_empty_structure();
obj.plot_properties = set_epsi_plot_properties;

% Add directory where .mat files exist
obj.Meta_Data.MATpath = pwd;

% Check if you're in a directory where .mat files exist
if ~exist(fullfile(obj.Meta_Data.MATpath,'Epsi_MATfile_TimeIndex.mat'),'file')
    error('Can''t load Epsi_MATfile_TimeIndex. Are you inside a directory with .mat files?')
end
%%
% Optional - add strings for mission, vehicle name, and deployment
obj.Meta_Data.mission = '';
obj.Meta_Data.vehicle_name = '';
obj.Meta_Data.deployment = '';

% Load the most recent .mat file in this directory and make a structure
% that looks like an epsi_class
tMaxOrig = now - 1;%%datenum(2021,7,3,0,0,0);
switch instrument
    case 'fctd'
        [obj,tMax] = getUpdatedData2(obj,tMaxOrig);
    case 'epsi'
        [obj,tMax] = getUpdatedData2(obj,tMaxOrig);
end


% Create an axes that will be the input for the first call to
% epsiPlot_epsi_ctd_alt_timeseries. All subsequent calls will reuse the set
% of axes created by that function.
ax = axes;
ax = epsiPlot_epsi_ctd_alt_timeseries(obj,0,ax);

% Initialize the timer to auto-update the realtime plot
EpsiPlot_time = timer;
EpsiPlot_time_Stop = false;

EpsiPlot_time.StartFcn = 'disp(''plotting epsi data...'');';
switch instrument
    case 'fctd'
        EpsiPlot_time.TimerFcn = [...
            'if EpsiPlot_time_Stop, '...
            'stop(EpsiPlot_time); '...
            'delete(EpsiPlot_time); '...
            'else, '...
            'disp([datestr(now)]); '...
            'try, '...
            '[obj,tMax] = getUpdatedData2(obj,tMax); '...
            'ax = fctdPlot_epsi_ctd_alt_timeseries(obj,0,ax); '...
            'catch err, '...
            'disp(err); '...
            'end; '...
            'end;'];
    case 'epsi'
        EpsiPlot_time.TimerFcn = [...
            'if EpsiPlot_time_Stop, '...
            'stop(EpsiPlot_time); '...
            'delete(EpsiPlot_time); '...
            'else, '...
            'disp([datestr(now)]); '...
            'try, '...
            '[obj,tMax] = getUpdatedData2(obj,tMax); '...
            'ax = epsiPlot_epsi_ctd_alt_timeseries(obj,0,ax); '...
            'catch err, '...
            'disp(err); '...
            'end; '...
            'end;'];
end
EpsiPlot_time.Period = 10;
EpsiPlot_time.BusyMode = 'drop';
EpsiPlot_time.Name = 'EpsiPlot_time';
EpsiPlot_time.Tag = 'EpsiPlot_time';
EpsiPlot_time.StopFcn = 'disp([datestr(now) '': Stopped EpsiPlot_time'']);';
EpsiPlot_time.ExecutionMode = 'fixedSpacing';
% EpsiPlot_time.ExecutionMode = 'singleShot';
EpsiPlot_time.TasksToExecute = Inf;
EpsiPlot_time.ErrorFcn = 'disp(''%%%%%%%%%%%%% Error %%%%%%%%%%%%%'');';

start(EpsiPlot_time);


