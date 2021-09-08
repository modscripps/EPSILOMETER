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
load(fullfile(Meta_Data.MATpath,'TimeIndex'))

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
            load([Meta_Data.MATpath '/' TimeIndex.filenames{myFileIdx} '.mat']);
        catch
            error(['Can''t load ' TimeIndex.filenames{myFileIdx}])
        end
    elseif length(myFileIdx)>1
        % Load the first file and name the structures 'Out'
        try
            load([Meta_Data.MATpath '/' TimeIndex.filenames{myFileIdx(1)} '.mat']);
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
        % Load the rest of the files, merging as you go (there shouldn't be
        % more than 2, so this should be pretty fast)
        for iF=2:length(myFileIdx)
            try
                load([Meta_Data.MATpath '/' TimeIndex.filenames{myFileIdx(iF)} '.mat']);
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
        clear epsiOut ctdOut altOut
        
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
        for iField=1:numel(ctdFields)
            Timeseries.ctd.(ctdFields{iField}) = ctd.(ctdFields{iField})(inRange);
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
        for iField=1:numel(epsiFields)
            Timeseries.epsi.(epsiFields{iField}) = epsi.(epsiFields{iField})(inRange);
        end
    end
    
    %% Add alt to Timeseries
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
    
    %% Add Meta_Data
    Timeseries.Meta_Data = Meta_Data;
    
else
    Timeseries = [];
end
