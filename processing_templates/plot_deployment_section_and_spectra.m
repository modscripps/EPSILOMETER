ccc

while 1<2
%% Define directories and depth range
lab_dir = '/Volumes/FCTD_EPSI/Deployments/0804_d2_ef1_pc2/';
processing_dir = '/Users/Shared/EPSI_PROCESSING/0804_d2_ef1_pc2/';
fctd_ctd_dir = '/Users/Shared/EPSI_PROCESSING/0804_d2_ef1_pc2/FCTD_mat';
depth_range = 700:1450;
% -----------------------------

% Meta_Data Process file
md = ['/Volumes/FCTD Softwares used in BLT 2022/EPSILOMETER_FCTD/'...
    'Meta_Data_Process/Meta_Data_Process_blt_2022.txt'];

% Rsync from Deployment raw directory to EPSI_PROCESSING raw directory
lab_dir_raw = fullfile(lab_dir,'raw');
processing_dir_raw = fullfile(lab_dir,'raw');
if ~exist(processing_dir_raw,'dir')
    mkdir(processing_dir_raw);
end
%com = sprintf('/usr/bin/rsync -auv %s %s',lab_dir_raw,processing_dir_raw);
com = sprintf('/usr/bin/rsync -auv %s %s',lab_dir,processing_dir);
unix(com);

%% Start epsi_class and read data with 'calc_micro' and 'make_FCTD'
% This will calculate epsilon and chi in the timeseries of each .mat file
% and feed those variables into FCTD .mat files
ec = epsi_class(processing_dir,md);
%ec.f_readData('calc_micro','make_FCTD',fctd_mat_dir)

%% Process new profiles
% This step will calculate epsilon and chi all over again (so yeah, this
% script is slooooow) but this time in the Profile structures.
ec.f_processNewProfiles;

%% Grid the profiles
ec.f_gridProfiles(depth_range);

%% Load the gridded profiles
load(fullfile(ec.Meta_Data.paths.profiles,'griddedProfiles'))

%% Plot temperature, epsilon, and chi with density contours
grid.bottom_depth=filloutliers(grid.bottom_depth,'linear');
close all
figure('units','inches','position',[0         0   15.3194   13.1111])
ax(1) = subtightplot(3,1,1);

pcolor(grid.dnum,grid.z,grid.t);
hold on
n_fill_bathy(grid.dnum,grid.bottom_depth)
[c,ch]=contour(grid.dnum,grid.z,grid.sgth,20,'k');
caxis([6 8])
shading flat
colorbar
n_fill_bathy(grid.dnum,grid.bottom_depth)
title('Temperature','fontname','Times New Roman','fontsize',20)
ylabel('Depth','fontname','Times New Roman','fontsize',20)

ax(2) = subtightplot(3,1,2);
pcolor(grid.dnum,grid.z,log10(grid.epsilon));
shading flat
colorbar
title('Epsilon')
ylabel('Depth','fontname','Times New Roman','fontsize',20);
hold on
[c,ch]=contour(grid.dnum,grid.z,grid.sgth,20,'k');
caxis([-10 -6])
n_fill_bathy(grid.dnum,grid.bottom_depth)

ax(3) = subtightplot(3,1,3);
pcolor(grid.dnum,grid.z,log10(grid.chi2));
shading flat
hold on
[c,ch]=contour(grid.dnum,grid.z,grid.sgth,20,'k');
n_fill_bathy(grid.dnum,grid.bottom_depth)
title('Chi 2')
colorbar
caxis([-10 -6])
ylabel('Depth','fontname','Times New Roman','fontsize',20)
xlabel(datestr(grid.dnum(1),"dd-mm"),'fontname','Times New Roman','fontsize',20)
[ax(:).YDir] = deal('reverse');
for a=1:length(ax)
    ax(a).XTick=grid.dnum(1):3/24:grid.dnum(end);
    ax(a).XTickLabel='';
end
ax(3).XTickLabel=datestr(grid.dnum(1):3/24:grid.dnum(end),'DD - HH:MM');
ax(3).XTickLabelRotation=45;

lp = linkprop([ax(:)],'xlim');

% Save section plot
save_name = fullfile(ec.Meta_Data.paths.figures,'deployment_sections');
eval(['export_fig ' save_name ' -png -r150 -nocrop']);
%print('-dpng2',save_name)


%% Plot 3 sets of spectra per profiles
%  - Loop through the profiles that haven't yet been plotted. Sort
%    epsilon_final by magnitude and plot spectra from the 10%, 50%, and
%    90% highest values
fig_list = dir(fullfile(ec.Meta_Data.paths.figures,'*spectra*.png'));
if ~isempty(fig_list)
    % The profile number is in characters 8-10 of the figure file name
    figNumCell = cellfun(@(C) C(8:10),{fig_list(:).name},'uniformoutput',0).';
else
    figNumCell = {''};
end

prof_list =  dir(fullfile(ec.Meta_Data.paths.profiles,'Profile*.mat'));
% The profile number is in characters #-# of the profile file name
profNumCell = cellfun(@(C) C(8:10),{prof_list(:).name},'uniformoutput',0).';

%  Find the profiles that don't yet have spectra plots
not_plotted = setdiff(profNumCell,figNumCell);

for iP=1:length(not_plotted)
    % Load the profile
    load(fullfile(ec.Meta_Data.paths.profiles,['Profile',not_plotted{iP}]));

    % Get the depths of the 10%, 50%, and 90% highest values of epsilon.
    [eps_sorted,idx_sorted] = sort(Profile.epsilon_final);
    eps_sorted_not_nan = eps_sorted(~isnan(eps_sorted));
    idx_sorted_not_nan = idx_sorted(~isnan(eps_sorted));
    n_eps = length(eps_sorted_not_nan);

    idx10 = idx_sorted_not_nan(round(n_eps*0.1));
    idx50 = idx_sorted_not_nan(round(n_eps*0.5));
    idx90 = idx_sorted_not_nan(round(n_eps*0.9));

    depths = [Profile.z(idx10);...
        Profile.z(idx50);...
        Profile.z(idx90)];

    % Loop through depths and make spectra  figures
    for iD=1:length(depths)

        plot_profile_and_spectra(Profile,depths(iD));

        % Save spectra plot
        save_name = fullfile(ec.Meta_Data.paths.figures,...
            sprintf('Profile%03.0f_%04.0fm_spectra',Profile.profNum,depths(iD)));
        eval(['export_fig ' save_name ' -png -r150 -nocrop']);
        %print('-dpng2',save_name)
        close all
    end
end


pause(60*13)
end %end while loop
