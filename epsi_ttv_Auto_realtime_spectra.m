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

% -------------------------------------------------------------------------
% --- USER CHOICES --------------------------------------------------------
totalSec = 120; %Number of seconds of data that will be stored in memory
plotSec = 30; %Number of seconds of data that will be plotted in timeseries
tscan = 4; %Length of scan in seconds
centerScan = tscan/2; %Plotted spectra will be centered tscan/2 seconds from the end of the timeseries
str_to_match = '*';
CTD_SN = '0237';

% Directory containing streaming raw data
rawDir = '/Users/Shared/FCTD_EPSI/RAW';
dirs.raw_incoming = rawDir;

% Choose time units
time_units = 'dnum'; %uncomment this if you set the datetime on SOM
%time_units = 'seconds'; %uncomment this if you did not set the datetime on SOM

% --- END USER CHOICES ----------------------------------------------------
% -------------------------------------------------------------------------


switch time_units
    case 'seconds'
        TMAX = 0;
    case 'dnum'
        TMAX = datenum(2023,11,1); %Will plot data after this starting point
end

% Read configuration data:
% First, try reading configuration data from the
% file. If that doesn't work, try reading from a
% % configuration file.
obj.Meta_Data.paths.raw_data = dirs.raw_incoming;
obj.Meta_Data.paths.data = dirs.raw_incoming; %Need a place to save Meta_Data (we don't need it, but epsiSetup_fill_meta_data saves it automatically)

% Is there a log csv file? Is there a config file? Or
% is config data inside the raw data files?
dir_has_log = dir(fullfile(dirs.raw_incoming,'Log*.csv'));
dir_has_config = dir(fullfile(dirs.raw_incoming,'*config*'));

spltpath=strsplit(path,':');
epsilib_path=spltpath{~cellfun(@isempty, ...
    cellfun(@(x) ...
    strfind(x,'epsilib'),spltpath, ...
    'UniformOutput',false))};
obj.Meta_Data.paths.process_library=fileparts(epsilib_path);

if ~isempty(dir_has_log) %if there is a log file...
    
    try
        obj.Meta_Data = create_metadata_from_deployment_log_v2(dir_has_log.name);
        obj.Meta_Data.AFE=obj.Meta_Data.epsi;
    catch err
        
        for j = 1:length(err.stack)
            disp([num2str(j) ' ' err.stack(j).name ' ' num2str(err.stack(j).line)]);
        end
        error('Failed to find config data (1)')
    end
    
elseif ~isempty(dir_has_config) %if there is a config file...
    
    try
        setup=mod_som_read_setup_from_config(dir_has_config.name);
    catch err
        for j = 1:length(err.stack)
            disp([num2str(j) ' ' err.stack(j).name ' ' num2str(err.stack(j).line)]);
        end
        error('Failed to find config data (2)')
        
    end
    % Fill Meta Data from setup data
    try
        obj.Meta_Data.CTD.SN = CTD_SN;
        setup.S49.sn = CTD_SN;
        obj.Meta_Data = epsiSetup_fill_meta_data(obj.Meta_Data,setup);
        
        fprintf('Meta_Data.paths.process_library is %s \n',obj.Meta_Data.paths.process_library);
        fprintf('Meta_Data.paths.data is %s \n',obj.Meta_Data.paths.data);
    catch err
        
        for j = 1:length(err.stack)
            disp([num2str(j) ' ' err.stack(j).name ' ' num2str(err.stack(j).line)]);
        end
        error('fill_meta_data failed (2)')
    end
    
else %if there is no log file or config file, look for config data inside the raw files
    % TODO 10/7/21 - Loop through more that just the
    % first file to look for $SOM3
    
    try
        setupfile=dir(fullfile(dirs.raw_incoming,'EPSI*'));
        setup=mod_som_read_setup_from_raw(fullfile(setupfile(1).folder,setupfile(1).name));
    catch err
        for j = 1:length(err.stack)
            disp([num2str(j) ' ' err.stack(j).name ' ' num2str(err.stack(j).line)]);
        end
        error(['Failed to read config data (3) - '...
            'this is often because mod_som_read_setup_from raw does not '...
            'have the correct offsets and lengths. When changes are made on '...
            'the hardware side, they have to be made here too.'])
    end
    % Fill Meta Data from setup data
    try
        obj.Meta_Data.CTD.SN = CTD_SN;
        setup.S49.sn = CTD_SN;
        obj.Meta_Data = epsiSetup_fill_meta_data(obj.Meta_Data,setup);
        
    catch err
        for j = 1:length(err.stack)
            disp([num2str(j) ' ' err.stack(j).name ' ' num2str(err.stack(j).line)]);
        end
        error('fill_meta_data failed (3)')
    end
    
end

%end

% Initialize obj with structures big enough to load at least one Epsi .mat
% file into (epsi, ctd, and alt strucutres)
obj.epsi = [];
obj = epsiSetup_make_empty_structure(obj);
obj.plot_properties = epsiSetup_set_plot_properties;
% Create Meta_Data
obj.Meta_Data.CTD.SN = CTD_SN;
setup.S49.sn = CTD_SN;
obj.Meta_Data = epsiSetup_fill_meta_data(obj.Meta_Data,setup);
obj.Meta_Data = epsiSetup_read_MetaProcess(obj.Meta_Data,...
    fullfile(obj.Meta_Data.paths.process_library,'Meta_Data_Process','Meta_Data_Process_blt_2022.txt'));

% Apply TMAX to structure tMax. Since the instruments sample at different
% rates, these will become slightly different from each other in the loop
% as new data come in.
field_list = {'epsi','ctd','alt','vnav','gps','ttv'};
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
    'try, '...
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
    '[~,~,~,ax] = ttv_plot_at_tMid(obj,tMid,tscan,plotSec,1,0,1,ax);, '...
    'catch err, '...
    'display_error_stack(err); '...
    'tStop = 1;'...
    'end; '...
    'end; '...
    'catch err, '...
    'display_error_stack(err); '...
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


