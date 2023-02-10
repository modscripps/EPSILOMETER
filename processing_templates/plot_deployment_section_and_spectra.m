ccc
while 1
%% Define some things
% -------------------------------------------
experiment_name     = 'TLC';
station_number      = 1;
deployment_number   = 1;
depl_dir            = '02_10_day4_epsi_d2';
depth_range         = 0:150;
minutes_between_processing = 3;
% -------------------------------------------

% Directories
%lab_dir = fullfile('/Volumes/FCTD_EPSI/Deployments/',depl_dir,'/');
processing_dir = fullfile('/Users/Shared/EPSI_PROCESSING/Processing',depl_dir);
%fctd_ctd_dir = '/Users/Shared/EPSI_PROCESSING/FCTD_mat';

% Meta_Data Process file
md = ['/Volumes/FCTD Softwares used in TLC 2023/EPSILOMETER_FCTD/'...
    'Meta_Data_Process/Meta_Data_Process_tlc_2023.txt'];
% % Rsync from Deployment raw directory to EPSI_PROCESSING raw directory
% lab_dir_raw = fullfile(lab_dir,'raw');
% processing_dir_raw = fullfile(processing_dir,'raw');
% if ~exist(processing_dir_raw,'dir')
%     mkdir(processing_dir_raw);
% end
% %com = sprintf('/usr/bin/rsync -auv %s %s',lab_dir_raw,processing_dir_raw);
% com = sprintf('/usr/bin/rsync -auv %s/ %s',lab_dir_raw,processing_dir_raw);
% unix(com);

%% Start epsi_class and read data with 'calc_micro' and 'make_FCTD'
% This will calculate epsilon and chi in the timeseries of each .mat file
% and feed those variables into FCTD .mat files
ec = epsi_class(processing_dir,md);

ec.Meta_Data.experiment = experiment_name;
ec.Meta_Data.deployment = deployment_number;
ec.Meta_Data.station = station_number;

%% Read all the data files and convert to mat
ec.f_readData

%% Process new profiles
% This step will calculate epsilon and chi all over again (so yeah, this
% script is slooooow) but this time in the Profile structures.
ec.f_processNewProfiles;

%% Grid the profiles
ec.f_gridProfiles(depth_range);

%% Load the gridded profiles
data=load(fullfile(ec.Meta_Data.paths.profiles,'griddedProfiles'));

%% Plot temperature, epsilon, and chi with density contours
data.GRID.bottom_depth=filloutliers(data.GRID.bottom_depth,'linear');
close all
figure('units','inches','position',[0         0   15.3194   13.1111])

dnummask=~isnan(data.GRID.dnum);
ax(1) =subtightplot(3,1,1);

pcolorjw(data.GRID.dnum(dnummask),data.GRID.z,data.GRID.t(:,dnummask));
hold on
%n_fill_bathy(data.GRID.dnum(dnummask),data.GRID.bottom_depth(dnummask))
%[c,ch]=contour(data.GRID.dnum(dnummask),data.GRID.z,data.GRID.sgth(:,dnummask),20,'k');
caxis([9 15])
cax1=colorbar;
grid(ax(1),'on');
%n_fill_bathy(data.GRID.dnum(dnummask),data.GRID.bottom_depth(dnummask))
title('Temperature','fontname','Times New Roman','fontsize',20)
ylabel('Depth','fontname','Times New Roman','fontsize',20)
ylabel(cax1,'Celsius','fontname','Times New Roman','fontsize',20)

ax(2) = subtightplot(3,1,2);
pcolorjw(data.GRID.dnum(dnummask),data.GRID.z,log10(data.GRID.epsilon_final(:,dnummask)));
cax2=colorbar;
title('Epsilon')
ylabel('Depth','fontname','Times New Roman','fontsize',20);
hold on
[c,ch]=contour(data.GRID.dnum(dnummask),data.GRID.z,data.GRID.sgth(:,dnummask),'k','levellist',-10:0.5:-6);
caxis([-10 -6])
grid(ax(2),'on');

%n_fill_bathy(data.GRID.dnum(dnummask),data.GRID.bottom_depth(dnummask))
ylabel(cax2,'W kg^{-1}','fontname','Times New Roman','fontsize',20)

ax(3) = subtightplot(3,1,3);
pcolorjw(data.GRID.dnum(dnummask),data.GRID.z,log10(data.GRID.chi1(:,dnummask)));
hold on
[c,ch]=contour(data.GRID.dnum(dnummask),data.GRID.z,data.GRID.sgth(:,dnummask),'k','levellist',-10:0.5:-5);
%n_fill_bathy(data.GRID.dnum(dnummask),data.GRID.bottom_depth(:,dnummask))
title('Chi 1')
cax3=colorbar;
grid(ax(3),'on');

caxis([-10 -5])
ylabel('Depth','fontname','Times New Roman','fontsize',20)
xlabel(datestr(data.GRID.dnum(1),"dd-mm"),'fontname','Times New Roman','fontsize',20)
[ax(:).YDir] = deal('reverse');
for a=1:length(ax)
    ax(a).XTick=data.GRID.dnum(1):2/24:data.GRID.dnum(end);
    ax(a).XTickLabel='';
end
ax(3).XTickLabel=datestr(data.GRID.dnum(1):2/24:data.GRID.dnum(end),'DD - HH:MM');
ax(3).XTickLabelRotation=45;
ylabel(cax3,'K^2 s^{-1}','fontname','Times New Roman','fontsize',20)

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
        save_name = fullfile(ec.Meta_Data.paths.figures,...
            sprintf('Profile%03.0f_%04.0fm_spectra',Profile.profNum,depths(iD)));
        eval(['export_fig ' save_name ' -png -r150 -nocrop']);
        %print('-dpng2',save_name)
        close all
    end %ALB end for ID 
    end %ALB end n_eps>0
end


pause(60*minutes_between_processing)
end %end while loop
