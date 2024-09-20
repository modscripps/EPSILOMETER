function [Timeseries] = epsiProcess_crop_timeseries(Meta_Data,tRange)
% Timeseries = crop_timeseries(Meta_Data,[tMin,tMax])
%
% Get a short piece of timeseries structure that you can use to compute
% turbulence variables.
%
% INPUTS:
%   tRange - range of [dnumMin dnumMax]
%   variableList (optional, to do) - default is all variables
%
% OUTPUT:
%   Timeseries - structure of epsi and ctd data to process turbulence
%   variables

%% Find the mat file(s) within the time range
load(fullfile(Meta_Data.paths.mat_data,'TimeIndex'))

% To do: First, decide if the input was in datenum or seconds since power up.

% Determine if tRange is given in seconds or dnum and define startFile and
% endFile accordingly
if tRange(1)>1e9 || tRange(1)<7e5
    tRangeChoice = 'time_s';
    startFile = find(tRange(1)>=TimeIndex.timeStart & tRange(1)<=TimeIndex.timeEnd);
    endFile = find(tRange(end)>=TimeIndex.timeStart & tRange(end)<=TimeIndex.timeEnd);
else
    tRangeChoice = 'dnum';
    startFile = find(tRange(1)>=TimeIndex.dnumStart & tRange(1)<=TimeIndex.dnumEnd);
    endFile = find(tRange(end)>=TimeIndex.dnumStart & tRange(end)<=TimeIndex.dnumEnd);
end
% If startFile is empty, start with the first file. If endFile is empty,
% end with the last file.
if isempty(startFile)
    startFile = 1;
end
if isempty(endFile)
    endFile = length(TimeIndex.dnumEnd);
end
myFileIdx = startFile:endFile;


if ~isempty(myFileIdx)
    % If the profile is contained in one file, load it. If more than one, load
    % and merge them.
    if length(myFileIdx)==1
        try
            MatData=load([Meta_Data.paths.mat_data '/' TimeIndex.filenames{myFileIdx} '.mat']);
            use MatData;
        catch
            error(['Can''t load ' TimeIndex.filenames{myFileIdx}])
        end
        use MatData;

    elseif length(myFileIdx)>1
        % Load the first file and name the structures 'Out'
        try
            MatData=load([Meta_Data.paths.mat_data '/' TimeIndex.filenames{myFileIdx(1)} '.mat']);
        catch
            error(['Can''t load ' TimeIndex.filenames{myFileIdx(1)}])
        end
        if isfield(MatData,'epsi')
            epsiOut = MatData.epsi;
            clear epsi
        end
        if isfield(MatData,'ctd')
            ctdOut = MatData.ctd;
            clear ctd
        end
        if isfield(MatData,'alt')
            altOut = MatData.alt;
            clear alt
        end
        if isfield(MatData,'isap')
            isapOut = MatData.isap;
            clear alt
        end
        if isfield(MatData,'vnav')
            vnavOut = MatData.vnav;
            clear vnav
        end
        if isfield(MatData,'ttv')
            ttvOut = MatData.ttv;
            clear ttv
        end
        if isfield(MatData,'gps')
            gpsOut = MatData.gps;
            clear ttv
        end
        % Load the rest of the files, merging as you go (there shouldn't be
        % more than 2, so this should be pretty fast)
        for iF=2:length(myFileIdx)
            try
                MatData=load([Meta_Data.paths.mat_data '/' TimeIndex.filenames{myFileIdx(iF)} '.mat']);
            catch
                error(['Can''t load ' TimeIndex.filenames{myFileIdx(iF)}])
            end
            if isfield(MatData,'epsi') && isstruct(MatData.epsi)
                epsiOut = epsiProcess_merge_mat_files(epsiOut,MatData.epsi);
            end
            if isfield(MatData,'ctd') && isstruct(MatData.ctd)
                ctdOut = epsiProcess_merge_mat_files(ctdOut,MatData.ctd);
            end
            if isfield(MatData,'alt') && isstruct(MatData.alt)
                altOut = epsiProcess_merge_mat_files(altOut,MatData.alt);
            end
            if isfield(MatData,'isap') && isstruct(MatData.isap)
                isapOut = epsiProcess_merge_mat_files(altOut,MatData.isap);
            end
            if isfield(MatData,'vnav') && isstruct(MatData.vnav)
                vnavOut = epsiProcess_merge_mat_files(vnavOut,MatData.vnav);
            end
            if isfield(MatData,'ttv') && isstruct(MatData.ttv)
                ttvOut = epsiProcess_merge_mat_files(ttvOut,MatData.ttv);
            end
            if isfield(MatData,'gps') && isstruct(MatData.gps)
                gpsOut = epsiProcess_merge_mat_files(gpsOut,MatData.gps);
            end
        end
        
        % Rename everything
        if isfield(MatData,'epsi')
            epsi = epsiOut;
        end
        if isfield(MatData,'ctd')
            ctd = ctdOut;
        end
        if isfield(MatData,'alt')
            alt = altOut;
        end
        if isfield(MatData,'isap')
            isap = isapOut;
        end
        if isfield(MatData,'vnav')
            vnav = vnavOut;
        end
        if isfield(MatData,'ttv')
            ttv = ttvOut;
        end
        if isfield(MatData,'gps')
            gps = gpsOut;
        end
        clear epsiOut ctdOut altOut vnavOut ttvOut isapOut gpsOut
        
    end
    
    %% Add filenames to Timeseries
    Timeseries.filenames = TimeIndex.filenames(myFileIdx);
    
    %% Add ctd to Timeseries
    if isfield(ctd,'dnum') || isfield(ctd,'time_s')
        switch tRangeChoice
            case 'dnum'
            inRange = ctd.dnum>=tRange(1) & ctd.dnum<=tRange(end);
            case 'time_s'
            inRange = ctd.time_s>=tRange(1) & ctd.time_s<=tRange(end);
        end
        
        ctdFields = fields(MatData.ctd);
        % Don't add any of the '_raw' fields
        notRaw = cell2mat(cellfun(@(C) isempty(strfind(C,'_raw')),ctdFields,'UniformOutput',0));
        ctdFields = ctdFields(notRaw);
        for iField=1:numel(ctdFields)
            if ~all(isnan(MatData.ctd.(ctdFields{iField})))  %NC 2/25/22 - Changed from checking for no nans to checking that it is not ALL nans
                Timeseries.ctd.(ctdFields{iField}) = ctd.(ctdFields{iField})(inRange);
            else
                Timeseries.ctd.(ctdFields{iField})=[];
            end
        end
    end
    
    %% Add epsi to Timeseries
    if isfield(epsi,'dnum') || isfield(epsi,'time_s')
        switch tRangeChoice
            case 'dnum'
            inRange = epsi.dnum>=tRange(1) & epsi.dnum<=tRange(end);
            case 'time_s'
            inRange = epsi.time_s>=tRange(1) & epsi.time_s<=tRange(end);
        end
        
        epsiFields = fields(epsi);
        % Don't add any of the '_count' fields
        notCount = cell2mat(cellfun(@(C) isempty(strfind(C,'_count')),epsiFields,'UniformOutput',0));
        epsiFields = epsiFields(notCount);
        for iField=1:numel(epsiFields)
            Timeseries.epsi.(epsiFields{iField}) = epsi.(epsiFields{iField})(inRange);
        end
    end
    
    %% Add alt to Timeseries
    if exist('alt','var')
    if isfield(alt,'dnum') || isfield(alt,'time_s')
        switch tRangeChoice
            case 'dnum'
            inRange = alt.dnum>=tRange(1) & alt.dnum<=tRange(end);
            case 'time_s'
            inRange = alt.time_s>=tRange(1) & alt.time_s<=tRange(end);
        end
        altFields = fields(alt);
        for iField=1:numel(altFields)
            Timeseries.alt.(altFields{iField}) = alt.(altFields{iField})(inRange);
        end
    end
    end
    %% Add isap to Timeseries
    if exist('isap','var')
    if isfield(isap,'dnum') || isfield(isap,'time_s')
        switch tRangeChoice
            case 'dnum'
            inRange = isap.dnum>=tRange(1) & isap.dnum<=tRange(end);
            case 'time_s'
            inRange = isap.time_s>=tRange(1) & isap.time_s<=tRange(end);
        end
        isapFields = fields(isap);
        for iField=1:numel(isapFields)
            Timeseries.isap.(isapFields{iField}) = isap.(isapFields{iField})(inRange);
        end
    end
    end


    %% Add vnav to Timeseries
    if exist('vnav','var')
        
        if isfield(vnav,'dnum') || isfield(vnav,'time_s')
            switch tRangeChoice
                case 'dnum'
                    inRange = vnav.dnum>=tRange(1) & vnav.dnum<=tRange(end);
                case 'time_s'
                    inRange = vnav.time_s>=tRange(1) & vnav.time_s<=tRange(end);
            end
            vnavFields = fields(vnav);
            for iField=1:numel(vnavFields)
                Timeseries.vnav.(vnavFields{iField}) = vnav.(vnavFields{iField})(inRange,:);
            end
        end
    end
    
    %% Add gps to Timeseries
    if exist('gps','var')
        
        if isfield(gps,'dnum') || isfield(gps,'time_s')
            switch tRangeChoice
                case 'dnum'
                    inRange = gps.dnum>=tRange(1) & gps.dnum<=tRange(end);
                case 'time_s'
                    inRange = gps.time_s>=tRange(1) & gps.time_s<=tRange(end);
            end
            gpsFields = fields(gps);
            for iField=1:numel(gpsFields)
                Timeseries.gps.(gpsFields{iField}) = gps.(gpsFields{iField})(inRange,:);
            end
        end
    end
    %% Add ttv to Timeseries
    if exist('ttv','var')
        if isfield(ttv,'dnum') || isfield(ttv,'time_s')
            switch tRangeChoice
                case 'dnum'
                    inRange = ttv.dnum>=tRange(1) & ttv.dnum<=tRange(end);
                case 'time_s'
                    inRange = ttv.time_s>=tRange(1) & ttv.time_s<=tRange(end);
            end
            ttvFields = fields(ttv);
            for iField=1:numel(ttvFields)
                if ~isstruct(ttv.(ttvFields{iField}))
                    Timeseries.ttv.(ttvFields{iField}) = ttv.(ttvFields{iField})(inRange,:);
                end
            end
        end
    end

    
    %% Add Meta_Data
    % ALB I am using the Meta_Data from the modmat.
    % I think it is fime to use any of them (many modmat per profiles) for a profile 
    try
        Timeseries.Meta_Data = MatData.Meta_Data;
    catch
        Timeseries.Meta_Data = [];
    end
    
else
    Timeseries = [];
end
