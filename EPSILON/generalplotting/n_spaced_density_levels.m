function [dens_levels] = n_spaced_density_levels(dens_mean,depth_mean,n_levels)

round_decimal = 1;

% Get rid of nans in mean density profile and then choose levels equally
% spaced in depth
nans = isnan(dens_mean);
dens_mean = dens_mean(~nans);
depth_mean = depth_mean(~nans);
idx = round(linspace(1,length(depth_mean),n_levels));
dens_levels = round(dens_mean(idx),round_decimal);