
% This is a wrapper script to set up input values and run
% epsiAuto_process_data.m
% 
% CREATE:
% input_struct, a structure containing the following fields
%
% Required:
%  .raw_dir = path to directory where the data are streaming in 
%  .process_dir = path to directory where data will be copied and where 
%              subdirectories raw, mat, and FCTDmat will be created
%  .Meta_Data_process_file = path to Meta_Data_Process text file
%
% Optional:
%  .str_to_match = (default '*'), String to search for files in raw_dir to copy. Often, we'll 
%                have data from the entire cruise and divide deployments up by 
%                day. You could select to only copy data from Feb 16th by 
%                setting this field to '*EPSI_B_PC2_22_02_16*'; 
%  .refresh_time_sec = (default 5*60), refresh period in seconds 
%  .version       = (default 4), version of mod_som_read_epsi_files.m to use 
% -------------------------------------------------------------------------

% Change these for each deployment
instrument = 'fctd';
input_struct.process_dir = '/Volumes/My Drive/DATA/TFO/RR2410/fctd_epsi/20240811_d7_CrossGS_SIOS_SION/';
input_struct.str_to_match = 'EPSI24';

% -------------------------------------------------------------------------

% These will probably be the same for the whole cruise
input_struct.raw_dir = '/Volumes/My Drive/DATA/TFO/RR2410/fctd_epsi/20240811_d7_CrossGS_SIOS_SION/raw/';
input_struct.Meta_Data_process_file = '~/ARNAUD/SCRIPPS/EPSILOMETER/Meta_Data_Process/MDP_tfo_2024.txt';
input_struct.refresh_time_sec =  5*60;
input_struct.cruise_specifics = 'tfo_2024';
switch instrument
    case 'epsi'
        input_struct.depth_array = 0:500;
    case 'fctd'
        input_struct.depth_array = 0:3000;
end

% fctdAuto_convert_raw_to_mat_and_plot(input_struct);
%
% automation to convert .raw files to .mat in epsi and fctd formats
% required by plotting function
%
% Nicole Couto adapted from autorunFastCTDConvert.m
% May 2021
% Updated: July 2022
%          May 2023 - removed plotting, this is now just an auto-processing
%          function
% -------------------------------------------------------------------------
% BEST PRACTICES:
%   Run this as a script inside a parent script that specifies options for
%   the data you wish to process (RUN_epsiAuto_process_data.m)
%
% INPUTS:
% input_struct contains the following fields
%
% Required:
%  .raw_dir = path to directory where the data are streaming in 
%  .process_dir = path to directory where data will be copied and where 
%              subdirectories raw, mat, and FCTDmat will be created
%  .Meta_Data_process_file = path to Meta_Data_Process text file
%
% Optional:
%  .str_to_match = (default '*'), String to search for files in raw_dir to copy. Often, we'll 
%                have data from the entire cruise and divide deployments up by 
%                day. You could select to only copy data from Feb 16th by 
%                setting this field to '*EPSI_B_PC2_22_02_16*'; 
%  .refresh_time_sec = (default 5*60), refresh period in seconds 
%  .version       = (default 4), version of mod_som_read_epsi_files.m to use 
% -------------------------------------------------------------------------
% Check inputs and apply defaults
field_list = fields(input_struct);
% Check that all required fields are defined
if ~any(contains(field_list,'raw_dir')) || ~any(contains(field_list,'process_dir')) || ~any(contains(field_list,'Meta_Data_process_file'))
  error('You must define ''raw_dir,'' ''process_dir,'' and ''Meta_Data_process_file'' as fields in the input structure')
else
  raw_dir = input_struct.raw_dir;
  process_dir = input_struct.process_dir;
  Meta_Data_process_file = input_struct.Meta_Data_process_file;
end
% Check for optional inputs and define them
if ~contains(field_list,'str_to_match')
  str_to_match = '*';
else
  str_to_match = input_struct.str_to_match;
end
if ~contains(field_list,'refresh_time_sec')
  refresh_time_sec = 5*60;
else
  refresh_time_sec = input_struct.refresh_time_sec;
end
if ~contains(field_list,'version ')
  version  = 4;
else
  version  = input_struct.version;
end
if ~contains(field_list,'starting_dnum')
    starting_dnum = now-14;
else
    starting_dnum = input_struct.starting_dnum;
end
if ~contains(field_list,'cruise_specifics')
    cruise_specifics = 'standard';
else
    cruise_specifics = input_struct.cruise_specifics;
end

% Define the directories for epsiProcess_convert_new_raw_to_mat
dirs.raw_incoming = raw_dir;
dirs.raw_copy  = fullfile(process_dir,'raw');
dirs.mat       = fullfile(process_dir,'mat');
experiment_dir = fileparts(process_dir);
dirs.fctd_cruise  = fullfile(experiment_dir,'FCTDmat');
dirs.fctd_deployment = fullfile(process_dir,'fctd_mat');
dirs.fctd_rot  = fullfile(experiment_dir,'FCTDrot');

% Create directories if they don't exist
if ~exist(process_dir,'dir')
    eval([ '!mkdir ' strrep(process_dir,' ','\ ')]);
end
if ~exist(dirs.raw_copy,'dir')
    eval([ '!mkdir ' strrep(dirs.raw_copy,' ','\ ')]);
end
if ~exist(dirs.mat,'dir')
    eval([ '!mkdir ' strrep(dirs.mat,' ','\ ')]);
end
if ~exist(dirs.fctd_cruise,'dir')
    eval([ '!mkdir ' strrep(dirs.fctd_cruise,' ','\ ')]);
end
if ~exist(dirs.fctd_deployment,'dir')
    eval([ '!mkdir ' strrep(dirs.fctd_deployment,' ','\ ')]);
end
if ~exist(dirs.fctd_rot,'dir')
    eval([ '!mkdir ' strrep(dirs.fctd_rot,' ','\ ')]);
end


% Copy the first file that matches str_to_match from raw_incoming into
% raw_copy - you need to have one file there for epsi_class to read the
% configuration information
if isfield(input_struct,'str_to_match')
    file_list = dir(fullfile(dirs.raw_incoming,[input_struct.str_to_match '*']));
else
    file_list = dir(fullfile(dirs.raw_incoming,'EPSI*'));
end
% eval(['!cp ' fullfile(file_list(1).folder,file_list(1).name) ' ' dirs.raw_copy]);

[SUCCESS,MESSAGE,MESSAGEID] = copyfile(fullfile(file_list(1).folder,file_list(1).name),dirs.raw_copy);
% % Initialize epsi_class in process_dir and create blank structures to fill
% % with data
% obj = epsi_class(process_dir,Meta_Data_process_file);
% obj = epsiSetup_make_empty_structure(obj);
% 
% field_list = {'epsi','ctd','alt','vnav','gps'};
% for iField=1:length(field_list)
%     tMax.(field_list{iField}) = starting_dnum;
% end
% % Create an axes that will be the input for the first call to
% % epsiPlot_epsi_ctd_alt_timeseries. All subsequent calls will reuse the set
% % of axes created by that function.
% ax = axes;
% ax = epsiPlot_timeseries(obj,0,ax);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ec = epsi_class(process_dir,Meta_Data_process_file);

disp([datestr(now) ': Converting new raw files to mat...']);
epsiProcess_convert_new_raw_to_mat(dirs,ec.Meta_Data,"noSync","doFCTD","fileStr",str_to_match,"cruise_specifics",cruise_specifics);
addpath(genpath('~/ARNAUD/SCRIPPS/FastCTD_MATLAB/'))
plot_fctd_sections;

