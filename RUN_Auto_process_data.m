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
<<<<<<< Updated upstream
=======
    % instrument = 'epsi';
>>>>>>> Stashed changes

end

newSurvey_flag=contains(str,'CTD.survey:');
if newSetup_flag
    surveyflag_str      = str(strfind(str,'CTD.survey:')+(0:100));
    surveyflag_str      = surveyflag_str(1:find(uint8(surveyflag_str)==10,1,'first'));
    surveyflag_name     = strsplit(surveyflag_str,':');
    surveylag_name      = surveyflag_name{2}(1:end-1);
    survey_name = surveylag_name;

else
    survey_name = input('What is the survey name? [d#_NAME]');

end

% ALB this routine now grab the survey name from the setup file 
% create the new process_dir
input_struct.process_dir = fullfile( ...
                           '/Users/Shared/EPSI_PROCESSING/TFO2024/Processed/', ...
                           survey_name); %This will create a directory with this name
if ~exist(fullfile(input_struct.process_dir,'raw'),'dir')
    mkdir(input_struct.process_dir);
end



% DO NOT CHANGE raw_dir 
input_struct.raw_dir = '/Users/Shared/EPSI_PROCESSING/TFO2024/Realtime_RAW/raw/';


% ALB getting rid of str_to_match
% looking for the folder "raw" in input_struct.process_dir
% and the lastest file in "raw" to get str_
if exist(fullfile(input_struct.process_dir,'raw'),'dir')
    listfile_process_dir=dir(fullfile(input_struct.process_dir,'raw'));
    if ~isempty(listfile_process_dir)
        input_struct.str_to_match = listfile_process_dir(end).name;
    else
        listfile_raw_dir=dir(input_struct.raw_dir);
        count=0;
        % quick look in the previous files if survey name exist.
        % In case blue matlab was launched few minutes after fctd_epsi
        while count<10
            path2setup1=fullfile(...
                listfile_raw_dir((length(listfile_raw_dir)-count)).folder,...
                listfile_raw_dir((length(listfile_raw_dir)-count)).name);
            fid=fopen(path2setup1,'r');
            fseek(fid,0,1);
            frewind(fid);
            str = fread(fid,'*char')';
            fclose(fid);
            surveyflag1_str      = str(strfind(str,'CTD.survey:')+(0:100));
            surveyflag1_str      = surveyflag1_str(1:find(uint8(surveyflag1_str)==10,1,'first'));
            surveyflag1_name     = strsplit(surveyflag1_str,':');
            surveylag1_name      = surveyflag1_name{2}(1:end-1);
            survey1_name         = surveylag1_name;

            if contains(survey1_name,survey_name)
                input_struct.str_to_match = ...
                    listfile_raw_dir(length(listfile_raw_dir)-count).name;
                count=count+1;
            else
                input_struct.str_to_match = ...
                    listfile_raw_dir(length(listfile_raw_dir)-count+1).name;
                % setting count to 10 to break out from while loop
                count=10; % The previous file is the first of the survey
            end % if survey1=survey
        end % end of while loop.
    end
else
    listfile_raw_dir=dir(input_struct.raw_dir);
    count=0;
    % quick look in the previous files if survey name exist.
    % In case blue matlab was launched few minutes after fctd_epsi
    while count<10
        path2setup1=fullfile(...
            listfile_raw_dir((length(listfile_raw_dir)-count)).folder,...
            listfile_raw_dir((length(listfile_raw_dir)-count)).name);
        fid=fopen(path2setup1,'r');
        fseek(fid,0,1);
        frewind(fid);
        str = fread(fid,'*char')';
        fclose(fid);
        surveyflag1_str      = str(strfind(str,'CTD.survey:')+(0:100));
        surveyflag1_str      = surveyflag1_str(1:find(uint8(surveyflag1_str)==10,1,'first'));
        surveyflag1_name     = strsplit(surveyflag1_str,':');
        surveylag1_name      = surveyflag1_name{2}(1:end-1);
        survey1_name         = surveylag1_name;

        if contains(survey1_name,survey_name)
            input_struct.str_to_match = ...
                listfile_raw_dir(length(listfile_raw_dir)-count).name;
            count=count+1;
        else
            input_struct.str_to_match = ...
                listfile_raw_dir(length(listfile_raw_dir)-count+1).name;
            % setting count to 10 to break out from while loop
            count=10; % The previous file is the first of the survey
        end % if survey1=survey
    end % end of while loop.
end


%ALB include survey name in the Meta_Data

% -------------------------------------------------------------------------

% These will probably be the same for the whole cruise
input_struct.Meta_Data_process_file = '/Volumes/Software_TFO2024/EPSILOMETER/Meta_Data_Process/MDP_tfo_2024.txt';
input_struct.refresh_time_sec =  2*60;
input_struct.cruise_specifics = 'tfo_2024';
switch instrument
    case {'epsi','EPSI'}
        input_struct.depth_array = 0:1600;
    case {'fctd','FCTD'}
        input_struct.depth_array = 0:2100;
end

% Set command window color
set_window_color('cyan')
% Run the processing script on a timer
switch instrument
    case {'epsi','EPSI'}
        epsiAuto_process_data
    case {'fctd','FCTD'}
        fctdAuto_process_data
end

