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
            load([Meta_Data.paths.mat_data '/' TimeIndex.filenames{myFileIdx} '.mat']);
        catch
            error(['Can''t load ' TimeIndex.filenames{myFileIdx}])
        end
    elseif length(myFileIdx)>1
        % Load the first file and name the structures 'Out'
        try
            load([Meta_Data.paths.mat_data '/' TimeIndex.filenames{myFileIdx(1)} '.mat']);
        catch
            error(['Can''t load ' TimeIndex.filenames{myFileIdx(1)}])
        end
        if exist('epsi','var')
            epsiOut = epsi;
            clear epsi
        end
        if exist('ctd','var')
            ctdOut = ctd;
            clear ctd
        end
        if exist('alt','var')
            altOut = alt;
            clear alt
        end
        if exist('vnav','var')
            vnavOut = vnav;
            clear vnav
        end
        if exist('fluor','var')
            fluorOut = fluor;
            clear fluor
        end

        % Load the rest of the files, merging as you go (there shouldn't be
        % more than 2, so this should be pretty fast)
        for iF=2:length(myFileIdx)
            try
                load([Meta_Data.paths.mat_data '/' TimeIndex.filenames{myFileIdx(iF)} '.mat']);
            catch
                error(['Can''t load ' TimeIndex.filenames{myFileIdx(iF)}])
            end
            if exist('epsi','var') && isstruct(epsi)
                epsiOut = epsiProcess_merge_mat_files(epsiOut,epsi);
            end
            if exist('ctd','var') && isstruct(ctd)
                ctdOut = epsiProcess_merge_mat_files(ctdOut,ctd);
            end
            if exist('alt','var') && isstruct(alt)
                altOut = epsiProcess_merge_mat_files(altOut,alt);
            end
            if exist('vnav','var') && isstruct(vnav)
                vnavOut = epsiProcess_merge_mat_files(vnavOut,vnav);
            end
            if exist('fluor','var') && isstruct(vnav)
                fluorOut = epsiProcess_merge_mat_files(fluorOut,fluor);
            end
            if exist('ttv','var') && isstruct(ttv)
                ttvOut = epsiProcess_merge_mat_files(ttvOut,ttv);
            end
        end
        
        % Rename everything
        if exist('epsi','var')
            epsi = epsiOut;
        end
        if exist('ctd','var')
            ctd = ctdOut;
        end
        if exist('alt','var')
            alt = altOut;
        end
        if exist('vnav','var')
            vnav = vnavOut;
        end
        if exist('fluor','var')
            fluor = fluorOut;
        end
        if exist('ttv','var')
            ttv = ttvOut;
        end
        clear epsiOut ctdOut altOut vnavOut fluorOut ttvOut
        
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
        
        ctdFields = fields(ctd);
        % Don't add any of the '_raw' fields
        notRaw = cell2mat(cellfun(@(C) isempty(strfind(C,'_raw')),ctdFields,'UniformOutput',0));
        ctdFields = ctdFields(notRaw);
        for iField=1:numel(ctdFields)
            if ~all(isnan(ctd.(ctdFields{iField})))  %NC 2/25/22 - Changed from checking for no nans to checking that it is not ALL nans
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

        %% Add fluor to Timeseries
    if exist('fluor','var')
        
        if isfield(fluor,'dnum') || isfield(fluor,'time_s')
            switch tRangeChoice
                case 'dnum'
                    inRange = fluor.dnum>=tRange(1) & fluor.dnum<=tRange(end);
                case 'time_s'
                    inRange = fluor.time_s>=tRange(1) & fluor.time_s<=tRange(end);
            end
            fluorFields = fields(fluor);
            for iField=1:numel(fluorFields)
                Timeseries.fluor.(fluorFields{iField}) = fluor.(fluorFields{iField})(inRange,:);
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
                Timeseries.ttv.(ttvFields{iField}) = ttv.(ttvFields{iField})(inRange,:);
            end
        end
    end

    
    %% Add Meta_Data
    Timeseries.Meta_Data = Meta_Data;
    
else
    Timeseries = [];
end
