%% Plot 3 sets of spectra per profiles
%  - Loop through the profiles that haven't yet been plotted. Sort
%    epsilon_final by magnitude and plot spectra from the 10%, 50%, and
%    90% highest values
fig_list = dir(fullfile(obj.Meta_Data.paths.figures,'*spectra*.png'));
if ~isempty(fig_list)
    % The profile number is in characters 8-10 of the figure file name
    figNumCell = cellfun(@(C) C(8:10),{fig_list(:).name},'uniformoutput',0).';
else
    figNumCell = {''};
end

prof_list =  dir(fullfile(obj.Meta_Data.paths.profiles,'Profile*.mat'));
% The profile number is in characters #-# of the profile file name
profNumCell = cellfun(@(C) C(8:10),{prof_list(:).name},'uniformoutput',0).';

%  Find the profiles that don't yet have spectra plots
not_plotted = setdiff(profNumCell,figNumCell);

for iP=1:length(not_plotted)
    % Load the profile
    load(fullfile(obj.Meta_Data.paths.profiles,['Profile',not_plotted{iP}]));

    % Get the depths of the 10%, 50%, and 90% highest values of epsilon.
    [eps_sorted,idx_sorted] = sort(Profile.epsilon_final);
    eps_sorted_not_nan = eps_sorted(~isnan(eps_sorted));
    idx_sorted_not_nan = idx_sorted(~isnan(eps_sorted));
    n_eps = length(eps_sorted_not_nan);

    if n_eps>0
    % ALB It chocked for a very small number of eps_sorted_not_nan 
    % ALB switching round to ceil so we have at leat neps*0.1=1
    idx10 = idx_sorted_not_nan(ceil(n_eps*0.1));
    idx50 = idx_sorted_not_nan(ceil(n_eps*0.5));
    idx90 = idx_sorted_not_nan(ceil(n_eps*0.9));

    depths = [Profile.z(idx10);...
        Profile.z(idx50);...
        Profile.z(idx90)];

    % Loop through depths and make spectra  figures
    for iD=1:length(depths)

        plot_profile_and_spectra(Profile,depths(iD));

        % Save spectra plot
        save_name = fullfile(obj.Meta_Data.paths.figures,...
            sprintf('Profile%03.0f_%04.0fm_spectra',Profile.profNum,depths(iD)));
        eval(['export_fig ' save_name ' -png -r150 -nocrop']);
        %print('-dpng2',save_name)
        close all
    end %ALB end for ID 
    end %ALB end n_eps>0
end