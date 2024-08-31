% RUN_epsiAuto_process_data.m
%
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
path2setup='~/ARNAUD/SCRIPPS/EPSILOMETER/acq/fctd_epsi_acq/build/fctd_epsi/Build/Products/Debug/Setup';
fid=fopen(path2setup,'r');
fseek(fid,0,1);
frewind(fid);
str = fread(fid,'*char')';
fclose(fid);
newSetup_flag=contains(str,'CTD.fishflag=');
if newSetup_flag
    fishflag_str      = str(strfind(str,'CTD.fishflag=')+(0:100));
    fishflag_str      = fishflag_str(1:find(uint8(fishflag_str)==10,1,'first'));
    fishflag_name      = strsplit(fishflag_str,'=');
    fishflag_name      = fishflag_name{2}(2:end-2);
    instrument = fishflag_name;

else
    % instrument = 'fctd';
    % instrument = 'fctd_tridente';
    instrument = 'epsi';

end

% Change these for each deployment
input_struct.process_dir = '/Users/Shared/FCTD_EPSI/RAW';
input_struct.str_to_match = 'EPSI24';

% -------------------------------------------------------------------------

% These will probably be the same for the whole cruise
input_struct.raw_dir = '/Users/Shared/FCTD_EPSI/RAW/raw';
input_struct.Meta_Data_process_file = '/Users/aleboyer/ARNAUD/SCRIPPS/EPSILOMETER/Meta_Data_Process/MDP_tfo_2024.txt';
input_struct.refresh_time_sec =  1*60;
input_struct.cruise_specifics = 'tfo_2024';
switch instrument
    case 'epsi'
        input_struct.depth_array = 0:500;
    case 'fctd'
        input_struct.depth_array = 0:3000;
end

% Set command window color
% set_window_color('cyan')
% Run the processing script on a timer
switch instrument
    case 'epsi'
        epsiAuto_process_data
    case 'fctd'
        fctdAuto_process_data
end

