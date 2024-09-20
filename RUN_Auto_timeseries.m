% RUN_epsiAuto_timeseries.m
%
%  This is a wrapper script to set up input values and run
%  epsiAuto_timeseries.m
%
% CREATE:
% input_struct, a structure containing the following fields:
% 
% Required:
%  .raw_dir = path to directory where the data are streaming in
%  .Meta_Data_process_file = path to Meta_Data_Process text file
%
% Optional:
%  .refresh_time_sec = refresh period in seconds (default 5)
%  .sec_to_store  = number of seconds of data that will be stored in memory
%                  (default 120 seconds)
%  .sec_to_plot   = number of seconds of data that will be plotted in
%                   timeseries (default 60)
%  .version       = version of mod_som_read_epsi_files.m to use (default 4)
%  .starting_dnum = earliest datenum to plot (default 2 weeks before today)
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% --- USER CHOICES --------------------------------------------------------
% --- From ftcd_epsi/Setup   -------------------------------------------------------
%TODO give a folder path to Setup
path2setup='/Volumes/Software_TFO2024/New_mod_fish_lib/MOD_fish_lib/EPSILOMETER/acq/fctd_epsi_acq/build/fctd_epsi/Build/Products/Debug/Setup';
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
    instrument = input('what fish are we using? [epsi,fctd]');

end



% Also plot spectra?
include_spectra = 0;

% Meta_Data process file (make sure this file has the correct serial
% numbers for CTD, s1, s2, t1, t2. If you're running fctd, you can leave
% s1, s2, t1, t2 = '115')
Meta_Data_process_file = 'MDP_tfo_2024.txt';

% These probably will be the same for the whole cruise

% input_struct.raw_dir = '/Users/Shared/EPSI_PROCESSING/TFO2024/Realtime_RAW/';
% input_struct.Meta_Data_process_file = '/Volumes/Software_TFO2024/EPSILOMETER/Meta_Data_Process/MDP_tfo_2024.txt';

input_struct.raw_dir = '/Users/Shared/EPSI_PROCESSING/TFO2024/Realtime_RAW/';
Meta_Data_process_dir = '/Volumes/Software_TFO2024/EPSILOMETER/Meta_Data_Process/';
input_struct.Meta_Data_process_file = fullfile(Meta_Data_process_dir,Meta_Data_process_file);
input_struct.refresh_time_sec = 2;

% --- END USER CHOICES ----------------------------------------------------
% -------------------------------------------------------------------------

% Set command window color
set_window_color('yellow')
close all;
% Run the realtime plotting script on a timer
switch instrument
    case {'epsi','EPSI'}
        if ~include_spectra
            epsiAuto_timeseries
        elseif include_spectra
            epsiAuto_timeseries_spectra
        end
    case {'fctd','FCTD'}
        fctdAuto_timeseries
    case 'fctd_tridente'
        fctdAuto_timeseries_tridente
end