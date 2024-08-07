addpath(genpath('/Volumes/Versailles/Epsidrive/EPSIlOMETER'))

%% Set deployment name 
% (this is the name of the directory where you saved the raw data)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
deployment_name = '1126_epsi14_minnow1';
epsi_number = 1;
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

nov_day = str2num(deployment_name(3:4));
switch epsi_number
    case 1
        if nov_day<19
            MDP_file = 'MDP_norse2023_minnow1_beforeNov19.txt';
        else
            MDP_file = 'MDP_norse2023_minnow1.txt';
        end

    case 2
        if nov_day<19
            MDP_file = 'MDP_norse2023_minnow2_beforeNov19.txt';
        else
            MDP_file = 'MDP_norse2023_minnow2.txt';
        end

end


%% Set ranges
% Set depth range for gridding
depth_array = 1:200;

% Set axes limits for plotting
limits.depth = [0 200];
limits.lon = [-9 -3];
limits.lat = [70 72];
limits.temp = [0 5];
limits.sal = [32 33];
limits.epsilon = [-10 -6];
limits.chi = [-10 -4];

%% Set paths
% Set path to metadata file
norse_meta_data = fullfile('/Volumes/Versailles/Epsidrive/EPSILOMETER/Meta_Data_Process/',MDP_file);

% Set path to save figures
fig_dir = '/Volumes/2023007019/DATA/MOD_profiling/epsi/PLOTS';

% Set path to save data
sharedrive_dir = fullfile('/Volumes/2023007019/DATA/MOD_profiling/epsi/',deployment_name);

% Set path to raw data  directory (where you uploaded data)
raw_dir = fullfile('/Volumes/Versailles/Epsidrive/raw_data/',deployment_name);

% Set path to processing directory (this script will copy data to here and
% process the data in this directory)
process_dir = fullfile('/Volumes/Versailles/Epsidrive/processing/',deployment_name);

%% Load GPS data
try
gps_file =  '/Volumes/2023007019/DATA/TS_surface/TSG_latest.mat';
load(gps_file)
catch
end

%% Make a copy of the raw data in the processing directory.
mkdir(process_dir)
mkdir(fullfile(process_dir,'raw'))
com = sprintf('/usr/bin/rsync -auv %s %s',fullfile(raw_dir,'*.modraw'),fullfile(process_dir,'raw'));
unix(com);

%% Copy a bench_config into the processing directory
eval(['copyfile  /Volumes/Versailles/Epsidrive/EPSILOMETER/config_files/bench_config ' fullfile(process_dir)]);

%% Read and process data
ec = epsi_class(process_dir,norse_meta_data);

%%
ec.f_readData;
ec.f_processNewProfiles;

%% Grid data
ec.f_gridProfiles(depth_array);

%% Plot data
ec.f_plot_epsi_map_and_sections(limits,[1,2]);

% Add GPS data to map
fig = gcf;
axM = fig.Children(end);
axes(axM)

g = ec.f_loadGrid;

try
    lon = interp1(TSG.dnum,TSG.lon,g.dnum);
    lat = interp1(TSG.dnum,TSG.lat,g.dnum);

    hold on
    scatter(lon,lat,60,g.dnum,'filled','markeredgecolor','none')
    clim([nanmin(g.dnum),nanmax(g.dnum)])
catch
end

% Adjust title
titleStr = axM.Title.String;
new_title = {titleStr; strrep(deployment_name,'_','\_')};
axM.Title.String = new_title;

% %% Copy figures to ship share drive.
% % Save section plot
% save_name = fullfile(fig_dir,deployment_name);
% eval(['savefig ' save_name]);
% eval(['export_fig ' save_name ' -png -r150 -nocrop']);
% 
% %% Copy data to ship share drive
% mkdir sharedrive_data_dir
% eval(['copyfile ' fullfile(process_dir,'*') ' ' sharedrive_data_dir]);


