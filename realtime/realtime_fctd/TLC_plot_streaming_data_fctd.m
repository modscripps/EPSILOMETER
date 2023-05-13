% REQUIRED INPUTS
% raw_dir = where data are streaming in
% away_dir = where we copy files to see the timeseries in real time
% processing_dir = where we copy files to process them for microstructure
% profiles, etc.

%% MANDATORY INPUTS -------------------------------
% Raw directory (where data is streaming in)
input_struct.raw_dir = '/Volumes/FCTD_EPSI/TLC2023/';
%input_struct.raw_dir = 'mod_admin@192.168.1.158:/Users/Shared/FCTD_EPSI/RAW/';

% Directory to copy raw files - CHANGE THIS FOR NEW DEPLOYMENT
input_struct.away_dir = ['/Users/Shared/EPSI_PROCESSING/Realtime/'...
                        '02_02_pumptest/'];
                    
% Directory for processing data - CHANGE THIS FOR NEW DEPLOYMENT                   
input_struct.processing_dir = ['/Users/Shared/EPSI_PROCESSING/Processing/'...
                        '02_02_pumptest/'];

% Meta Data file
input_struct.Meta_Data_process_file = ...
    ['/Volumes/FCTD Softwares used in TLC 2023/EPSILOMETER_FCTD/'...
    'Meta_Data_Process/Meta_Data_Process_tlc_2023.txt'];
% FCTD directories 
input_struct.FCTDmat = '/Users/Shared/EPSI_PROCESSING/FCTD_mat/mat/';
input_struct.FCTDrot = '/Users/Shared/EPSI_PROCESSING/FCTD_mat/rot/';

%% OPTIONAL INPUTS -------------------------------
input_struct.str_to_match = 'EPSI23_02_02_2207'; %CHANGE THIS FOR NEW DEPLOYMENT. (Don't use *, this is going into a strfind)
input_struct.refresh_time_sec = 3;
%input_struct.version = 4;
%input_struct.starting_dnum = datenum(2023,1,1);



% 
% % REQUIRED INPUTS
% % raw_dir = where data are streaming in
% % away_dir = where we copy files to see the timeseries in real time
% % processing_dir = where we copy files to process them for microstructure
% % profiles, etc.
% input_struct.raw_dir = '/Users/Shared/EPSI_PROCESSING/RAW/';
% %input_struct.raw_dir = 'mod_admin@192.168.1.158:/Users/Shared/FCTD_EPSI/RAW/';
% input_struct.away_dir = ['/Users/Shared/EPSI_PROCESSING/Realtime/'...
%                         '0812_sta08_d17_ef1_pc2/']; %CHANGE THIS FOR NEW DEPLOYMENT
% input_struct.processing_dir = ['/Users/Shared/EPSI_PROCESSING/Processing/'...
%                         '0812_sta08_d17_ef1_pc2/']; %CHANGE THIS FOR NEW DEPLOYMENT
% input_struct.Meta_Data_process_file = ...
%     ['/Volumes/FCTD Softwares used in BLT 2022/EPSILOMETER_FCTD/'...
%     'Meta_Data_Process/Meta_Data_Process_blt_2022.txt'];
% 
% % OPTIONAL INPUTS
% input_struct.str_to_match = 'EPSI22_08_12_17'; %CHANGE THIS FOR NEW DEPLOYMENT. (Don't use *, this is going into a strfind)
% input_struct.refresh_time_sec = 5;
% %input_struct.version = 4;
% input_struct.starting_dnum = datenum(2022,7,31);
% 
% % FCTD directories
% input_struct.FCTDmat = '/Users/Shared/EPSI_PROCESSING/FCTD_mat/mat/';
% input_struct.FCTDrot = '/Users/Shared/EPSI_PROCESSING/FCTD_mat/rot/';
% 
fctdAuto_convert_raw_to_mat_and_plot
