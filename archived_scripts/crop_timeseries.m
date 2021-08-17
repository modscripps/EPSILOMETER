
function [Timeseries] = crop_timeseries(ec,Meta_Data,tRange)
% Timeseries = crop_timeseries(Meta_Data,[tMin,tMax])
%
% Get a short piece of timeseries structure that you can use to compute
% turbulence variables.
%
% INPUTS:
%   Meta_Data
%   tRange - range of [tMin tMax]
%
% OUTPUT:
%   Timeseries - structure of epsi and ctd data to process turbulence
%   variables

%% Load epsi and ctd
epsi = ec.epsi;
ctd = ec.ctd;

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
        case {'efe_badblocks','efe_time_sendout','efe_time_laptop', ...
              'efe_bad_blocks','efe_block_start'}
            
        otherwise
            Timeseries.(epsiFields{iField}) = epsi.(epsiFields{iField})(inRange);   
    end
end

%% Smooth the pressure channel and clean up
if isfield(Timeseries,'P')
    Timeseries.P = filloutliers(Timeseries.P,'center','movmedian',1000);
end
Timeseries = structfun(@(x) fillmissing(x,'linear'),Timeseries,'Un',0);

Timeseries.profNum = 0;