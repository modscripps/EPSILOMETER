function [] = epsiProcess_update_PressureTimeseries(Meta_Data,dirname,ctd,profile_dir)

if ~exist(fullfile(dirname, 'PressureTimeseries.mat'),'file')
    %Make a big array to fill CTD pressure and time. (24 hours)*(3600
    %sec/hour)*(16 Hz). If the deployment ends up being longer than 24
    %hours, the array will keep growing at every iteration and can slow
    %down processing speed.
    
    % ALB TODO use Meta_Data to define array length.
    arrayLength = 24*3600*16;

    % If it's the first time, make a big array of 'emptyFlag'
    PressureTimeseries.dnum = nan(arrayLength,1);
    PressureTimeseries.time_s = nan(arrayLength,1);
    PressureTimeseries.P = nan(arrayLength,1);
    idxLast = 0;
    nRecords = length(ctd.P);

    if isfield(PressureTimeseries,'dnum')
        PressureTimeseries.dnum(idxLast+1:idxLast+nRecords) = ctd.dnum;
    end
    PressureTimeseries.time_s(idxLast+1:idxLast+nRecords) = ctd.time_s;
    PressureTimeseries.P(idxLast+1:idxLast+nRecords) = ctd.P;

elseif exist(fullfile(dirname, 'PressureTimeseries.mat'),'file')
    % If PressureTimeseries already exists, find the new values of
    % dnum, and add them to the end of the record.
    load(fullfile(dirname, 'PressureTimeseries.mat'));
    idxLast = find(~isnan(PressureTimeseries.dnum),1,'last');
    maxValue = nanmax(PressureTimeseries.dnum(1:idxLast));

    newDnum = ctd.dnum(ctd.dnum>maxValue);
    newTime = ctd.time_s(ctd.dnum>maxValue);
    newP = ctd.P(ctd.dnum>maxValue);
    nRecords = length(newDnum);

    PressureTimeseries.dnum(idxLast+1:idxLast+nRecords) = newDnum;
    PressureTimeseries.time_s(idxLast+1:idxLast+nRecords) = newTime;
    PressureTimeseries.P(idxLast+1:idxLast+nRecords) = newP;

end

% Define (or redefine) the upcasts and downcasts
[PressureTimeseries] = epsiProcess_get_profiles_from_PressureTimeseries(PressureTimeseries,Meta_Data);
save(fullfile(dirname, 'PressureTimeseries.mat'),'PressureTimeseries');

end
