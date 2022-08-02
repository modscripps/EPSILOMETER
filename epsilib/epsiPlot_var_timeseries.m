function [var_timeseries] = epsiPlot_var_timeseries(obj,var_name)
% [var_timeseries] = epsiPlot_var_timeseries(obj,var_name)
% Plot the full timeseries of a single variable from an epsi_class object
%
% INPUTS
%   obj - epsi_class object
%   var_name - name of variable (e.g. t1_volt,s2_volt,T,dPdt)
%
% OUTPUTS
%   var_timeseries.dnum - datenum array
%   var_timeseries.data - data array
%

% Get the data
[var_timeseries] = epsiProcess_get_var_timeseries(obj,var_name);

% Plot the data
figure
plot(var_timeseries.dnum,var_timeseries.data,'.');
datetick('x','keeplimits')
title(strrep(var_name,'_','\_'))