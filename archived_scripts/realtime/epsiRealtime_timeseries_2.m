% epsiRealtime_timeseries_2.m - call from autorunEpsiConvert so you don't
% have to load .mat files
obj.epsi = epsi;
obj.ctd = ctd;
obj.alt = alt;
tMax = nanmax(data.epsi.epsidnum);


% Initialize obj with structures big enough to load at least one Epsi .mat
% file into (epsi, ctd, and alt strucutres)
%obj = make_empty_structure();
obj.plot_properties = set_epsi_plot_properties;

% Add directory where .mat files exist
obj.Meta_Data.MATpath = pwd;

% Optional - add strings for mission, vehicle name, and deployment
obj.Meta_Data.mission = '';
obj.Meta_Data.vehicle_name = '';
obj.Meta_Data.deployment = '';

% % Load the most recent .mat file in this directory and make a structure
% % that looks like an epsi_class
% tMaxOrig = datenum(2021,6,30,0,0,0);
% [obj,tMax] = getUpdatedData(obj,tMaxOrig);

% Create an axes that will be the input for the first call to
% epsiPlot_epsi_ctd_alt_timeseries. All subsequent calls will reuse the set
% of axes created by that function.
ax = axes;
ax = epsiPlot_epsi_ctd_alt_timeseries(obj,0,ax);

% Initialize the timer to auto-update the realtime plot
EpsiPlot_time = timer;
EpsiPlot_time_Stop = false;

EpsiPlot_time.StartFcn = 'disp(''plotting epsi data...'');';
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


