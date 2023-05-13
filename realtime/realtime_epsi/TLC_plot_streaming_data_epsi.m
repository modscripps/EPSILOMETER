% REQUIRED INPUTS
% raw_dir = where data are streaming in
% away_dir = where we copy files to see the timeseries in real time
% processing_dir = where we copy files to process them for microstructure
% profiles, etc.

%% MANDATORY INPUTS -------------------------------
% Raw directory (where data is streaming in)
%input_struct.raw_dir = '/Volumes/FCTD_EPSI/TLC2023/';
input_struct.raw_dir = '/Users/ncouto/Desktop/RAW';
%input_struct.raw_dir = 'mod_admin@192.168.1.158:/Users/Shared/FCTD_EPSI/RAW/';

% Directory to copy raw files - CHANGE THIS FOR NEW DEPLOYMENT
%input_struct.away_dir = ['/Users/Shared/EPSI_PROCESSING/Realtime/'...
%                        '02_10_day4_epsi_d2/'];
input_struct.away_dir = ['/Users/ncouto/Desktop/test_simulator/realtime/02_10_day4_epsi_d2/'];
                    
% Directory for processing data - CHANGE THIS FOR NEW DEPLOYMENT                   
%input_struct.processing_dir = ['/Users/Shared/EPSI_PROCESSING/Processing/'...
%                        '02_10_day4_epsi_d2/'];
input_struct.processing_dir = ['/Users/ncouto/Desktop/test_simulator/processing/02_10_day4_epsi_d2/'];

% Meta Data file
input_struct.Meta_Data_process_file = ...
    ['/Users/ncouto/GitHub/EPSILOMETER/'...
    'Meta_Data_Process/Meta_Data_Process_tlc_2023.txt'];
% FCTD directories 
input_struct.FCTDmat = '/Users/ncouto/Desktop/test_simulator/FCTD_mat/mat/';
input_struct.FCTDrot = '/Users/ncouto/Desktop/test_simulator/FCTD_mat/rot/';

%% OPTIONAL INPUTS -------------------------------
input_struct.str_to_match = 'EPSI23_02_10'; %CHANGE THIS FOR NEW DEPLOYMENT. (Don't use *, this is going into a strfind)
input_struct.refresh_time_sec = 3;
%input_struct.version = 4;
%input_struct.starting_dnum = datenum(2023,1,1);

%% Make directories for copying and processing raw data
if ~exist(fullfile(input_struct.away_dir),'dir')
    mkdir(fullfile(input_struct.away_dir,'raw'));
end
if ~exist(fullfile(input_struct.processing_dir),'dir')
    mkdir(fullfile(input_struct.processing_dir,'raw'));
end

%% Run the script
epsiAuto_convert_raw_to_mat_and_plot
