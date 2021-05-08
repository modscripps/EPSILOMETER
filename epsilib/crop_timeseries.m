function [Timeseries] = crop_timeseries(obj,tRange)
% Timeseries = crop_timeseries(obj,tRange)
%
% Get a short piece of timeseries structure that you can use to compute
% turbulence variables.
%
% INPUTS:
%   obj - epsi_class object
%   tRange - range of epsitime [tMin tMax]
%
% OUTPUT:
%   Timeseries - structure of epsi and ctd data to process turbulence
%   variables


%% Load ctd and epsi data
Meta_Data = obj.Meta_Data;

load(fullfile(Meta_Data.CTDpath,['ctd_' Meta_Data.deployment]));
load(fullfile(Meta_Data.Epsipath,['epsi_' Meta_Data.deployment]));

%% Add ctd to Timeseries
if isfield(ctd,'ctdtime')
    inRange = ctd.ctdtime>=tRange(1) & ctd.ctdtime<=tRange(end);
    
    ctdFields = fields(ctd);
    for iField=1:numel(ctdFields)
        switch ctdFields{iField}
            case 'timestamp'
                Timeseries.ctd_timestamp = ctd.(ctdFields{iField})(inRange);
            otherwise
                Timeseries.(ctdFields{iField}) = ctd.(ctdFields{iField})(inRange);
        end
    end
end

%% Add epsi to Timeseries
inRange = epsi.epsitime>=tRange(1) & epsi.epsitime<=tRange(end);

epsiFields = fields(epsi);
for iField=1:numel(epsiFields)
    switch epsiFields{iField}
        case 'timestamp'
            Timeseries.epsi_timestamp = epsi.(epsiFields{iField})(inRange);
        case {'efe_badblocks','efe_time_sendout','efe_time_laptop'}
            
        otherwise
            Timeseries.(epsiFields{iField}) = epsi.(epsiFields{iField})(inRange);   
    end
end

%% Smooth the pressure channel and clean up
if isfield(Timeseries,'P')
    Timeseries.P = filloutliers(Timeseries.P,'center','movmedian',1000);
end
Timeseries = structfun(@(x) fillmissing(x,'linear'),Timeseries,'Un',0);