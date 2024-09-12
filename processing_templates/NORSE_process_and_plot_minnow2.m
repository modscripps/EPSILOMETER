%% Set deployment name 
% (this is the name of the directory where you saved the raw data)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
deployment_name = '1109_lab_tests';

% Serial numbers
sn.SBE = '0237';
sn.s1 = 353;
sn.s2 = 351;
sn.t1 = 243;
sn.t2 = 274;
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% Set ranges
% Set depth range for gridding
depth_range = [0 200];

% Set axes limits for plotting
limits.depth = [0 200];
limits.lon = [-10 -2];
limits.lat = [69 72];
limits.temp = [0 6];
limits.sal = [30 34];
limits.epsilon = [-10 -6];
limits.chi = [-10 -6];

%% Set paths
% Set path to metadata file
norse_meta_data = '/Volumes/Epsidrive/EPSILOMETER/Meta_Data_Process/MDP_norse2023.txt';

% Set path to save figures
fig_dir = '/Volumes/2023007019/DATA/MOD_profiling/epsi/PLOTS';

% Set path to raw data  directory (where you uploaded data)
raw_dir = fullfile('/Volumes/Epsidrive/raw_data/',deployment_name);

% Set path to processing directory (this script will copy data to here and
% process the data in this directory)
process_dir = fullfile('/Volumes/Epsidrive/processing/',deployment_name);

%% Make a copy of the raw data in the processing directory.
com = sprintf('/usr/bin/rsync -auv %s %s',raw_dir,fullfile(process_dir,'raw'));
unix(com);

%% Copy a bench_config into the processing directory
copyfile /Volumes/Epsidrive/EPSILOMETER/config_files/bench_config fullfile(process_dir);


%% Read and process data
ec = epsi_class(process_dir,norse_meta_data);
ec.f_readData;
ec.f_processNewProfiles;
ec.f_gridProfiles(depth_range);

%% Plot data
ec.f_plot_epsi_map_and_sections(limits);

%% Copy figures to ship share drive.
% Save section plot
save_name = fullfile(fig_dir,'deployment_name');
eval(['savefig ' save_name ' -png -r150 -nocrop']);
eval(['export_fig ' save_name ' -png -r150 -nocrop']);


