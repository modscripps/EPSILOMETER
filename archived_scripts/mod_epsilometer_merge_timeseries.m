function [Timeseries]=mod_epsilometer_merge_timeseries(Meta_Data,tRange)

% [Timeseries]=mod_epsilometer_merge_timeseries(Meta_Data);
%
% Merge ctd and epsi timeseries into one structure


%% Load ctd and epsi data
load(fullfile(Meta_Data.CTDpath,['ctd_' Meta_Data.deployment]));
load(fullfile(Meta_Data.Epsipath,['epsi_' Meta_Data.deployment]));

%% Add ctd to Timeseries
ctdFields = fields(ctd);
for iField=1:numel(ctdFields)
    switch ctdFields{iField}
        case 'timestamp'
            Timeseries.ctd_timestamp = ctd.(ctdFields{iField});
        otherwise
            Timeseries.(ctdFields{iField}) = ctd.(ctdFields{iField});
    end
end

%% Add epsi to Timeseries
epsiFields = fields(epsi);
for iField=1:numel(epsiFields)
    switch epsiFields{iField}
        case 'timestamp'
            Timeseries.epsi_timestamp = epsi.(epsiFields{iField});
        otherwise
            Timeseries.(epsiFields{iField}) = epsi.(epsiFields{iField});
            
    end
end

%% Smooth the pressure channel and clean up
Timeseries.P = filloutliers(Timeseries.P,'center','movmedian',1000);
Timeseries = structfun(@(x) fillmissing(x,'linear'),Timeseries,'Un',0);
