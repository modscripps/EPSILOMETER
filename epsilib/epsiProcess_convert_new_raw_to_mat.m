function [matData] = epsiProcess_convert_new_raw_to_mat(dirs,Meta_Data,varargin)
% epsiProcess_convert_new_raw_to_mat
%   Converts new raw files from dirs.raw_incoming to mat files and saves them in dir.raw_copy
%
% aleboyer@ucsd.edu adding ttv and fluor
% Nicole Couto adapted from FCTD_MakeMatFromRaw.m
% May-July 2021
% April 2022
% --------------------------------------------------------
% INPUTS:
%   dirs.raw_incoming  - where raw data is streaming in
%       .raw_copy      - where raw data should be copied for the current
%                        deployment
%       .mat           - directory for .mat files converted from raw files
%       .fctd_cruise   - directory for .mat files formatted for FCTD
%                        processing, all files saved for whole cruise
%       .fctd_deployment - directory for .mat files formatted for FCTD
%                        processing, files just for current deployment
%   Meta_Data          - epsi Meta_Data structure
%
% OPTIONAL INPUTS:
%   'noSync'
%       - don't rsync files from RawDir to RawDir Away
%   'doFCTD'
%       - make FCTD files compatible with MOD FCTD processing library
%   'calc_micro'
%       - true/false to calculate microstructure data
%   'fileStr', 'string_to_match'
%       - use this flag along with a search string
%       if you want to limit rsyncing files to only those that match a certain
%       pattern. This is useful if you're storing all the raw data from a
%       cruise in RawDir, but you want to separate deployments into different
%       RawDirAway. You could specify only files from July 17th as
%       'fileStr','*07_17*com =

%    'version, version_number
%       - specify the version number (1,2,3,4) of mod_som_read_epsi_files.m
%    'display_file_data'
%       - display size of file, time, pressure, and altimeter
%    'cruise_specifics'
%       - on BLT 2021 cruise, we output microconductivity and fluorometer
%       data from epsi channels to the FCTD .mat structure
%
% Written by Jody Klymak
% Updated 2011 06 21 by San Nguyen
% Updated 2012 09 29 by San Nguyen for EquatorMix2012
% Adapted 2021-23 by Nicole Couto for Epsilometer data
% Adapted 2024-08-08 by aleboyer for ISA500

% NC - Make matData for output even if there is no new data
matData.epsi  = [];
matData.ctd   = [];
matData.alt   = [];
matData.isap  = [];
matData.act   = [];
matData.vnav  = [];
matData.gps   = [];
matData.seg   = [];
matData.spec  = [];
matData.micro = [];
matData.ttv   = [];
matData.fluor = [];

% NC - Only rsync files with the desired suffix
suffixStr = Meta_Data.rawfileSuffix; %ex. .raw, .ascii, etc
suffixSearch = ['*' suffixStr];

%% Set default parameters
% (If you specified them in call to this function)they will be changed in
% the next step
rSync      = true;  % Copy raw files from dirs.raw_incoming to dirs.raw_copy
doFCTD     = false; % Make separate .mat files in dirs.fctd* readable by FCTD processing scripts
calc_micro = false; % Calculate epsilon, chi, for timeseries (before dividing into profiles)
version    = 4;     % Default version of mod_som_read_epsi_files
display_file_data_flag    = false; % By default, DON'T display file information to the screen
cruise_specifics = 'standard'; % By default, DON'T add microconductivity and fluorometer data to FCTD structure
fileStr    = false;

%% Loop through the number of varargin arguments

argsNameToCheck = {'noSync',...             %1
    'doFCTD',...             %2
    'calc_micro',...         %3
    'fileStr',...            %4
    'version',...            %5
    'display_file_data',...  %6
    'cruise_specifics'};     %7

index = 1; %Initialize index of argsNameToCheck
% Number of items remaining (this is the number of argsNameToCheck minus
% the number of extra parameters that go with the arguments. For example,
% 'version' expects another parameter that follows it: 'version', 3.
% Similarly, 'fileStr' expects a string after it: 'fileStr',
% 'EPSI_22_04_12*'
n_items = nargin-2;

while (n_items > 0)
    argsMatch = strcmpi(varargin{index},argsNameToCheck);
    i = find(argsMatch,1);
    if isempty(i)
        error('MATLAB:epsiProcess_convert_new_raw_to_mat:wrongOption','Incorrect option specified: %s', varargin{index});
    end

    switch i
        case 1 % noSync
            rSync = false;
            index = index +1;
            n_items = n_items-1;
        case 2 %doFCTD
            idx_char=find(cellfun(@(x) (ischar(x)|isstring(x)),varargin));
            idxFlag = idx_char(cell2mat(cellfun(@(C) contains(C,'doFCTD'),varargin(idx_char),'uniformoutput',0)));
            doFCTD = true;
            index = index+1;
            n_items = n_items-1;
        case 3 %calc_micro
            % Find the index of varargin that = 'calc_micro'. The following
            % index contains the version number
            idx_char=find(cellfun(@(x) (ischar(x)|isstring(x)),varargin));
            idxFlag = idx_char(cell2mat(cellfun(@(C) contains(C,'calc_micro'),varargin(idx_char),'uniformoutput',0)));
            calc_micro = varargin{idxFlag+1};
            index = index+2;
            n_items = n_items-2;
        case 4 %fileStr % NC - added optional fileStr input
            fileStr = true;
            % Find the index of varargin that = 'fileStr'. The following
            % index contains 'str_to_match'
            idx_char=find(cellfun(@(x) (ischar(x)|isstring(x)),varargin));
            idxFlag = idx_char(cell2mat(cellfun(@(C) contains(C,'fileStr'),varargin(idx_char),'uniformoutput',0)));
            str_to_match = varargin{idxFlag+1};
            index = index+2; %+2 because the following varargin will str_to_match
            n_items = n_items-2;
        case 5 %version % NC 10/11/21 - added optional version input
            % Find the index of varargin that = 'version'. The following
            % index contains the version number
            idx_char=find(cellfun(@(x) (ischar(x)|isstring(x)),varargin));
            idxFlag = idx_char(cell2mat(cellfun(@(C) contains(C,'version'),varargin(idx_char),'uniformoutput',0)));
            version = varargin{idxFlag+1};
            index = index+2; %+2 because the following varargin will be the version number
            n_items = n_items-2;
        case 6 %display_file_data
            display_file_data_flag = true;
            index = index+1;
            n_items = n_items-1;
        case 7 %cruise_specifics
            idx_char=find(cellfun(@(x) (ischar(x)|isstring(x)),varargin));
            idxFlag = idx_char(cell2mat(cellfun(@(C) contains(C,'cruise_specifics'),varargin(idx_char),'uniformoutput',0)));
            cruise_specifics = varargin{idxFlag+1};
            index = index+2; %+2 because the following varargin will be the version number
            n_items = n_items-2;
    end
end

%% Define directories and check that they exist
% Incoming raw files
if ~exist(dirs.raw_incoming,'dir')
    error('Cannot find raw directory (dirs.raw_incoming): %s',RawDirDuplicate);
end
% Mat files, converted from raw
MatDir = dirs.mat;
if ~exist(MatDir,'dir')
    error('Cannot find .mat directory (dirs.mat): %s',MatDir);
end
% Copy of raw files
if rSync
    if ~exist(dirs.raw_copy,'dir')
        error('Cannot find duplicate raw directory (dirs.raw_copy): %s',RawDir);
    end
end
% Mat files, formatted for FCTD processing
if doFCTD
    if ~exist(dirs.fctd_cruise,'dir') % Bethan 20 June 2021: Added FCTD directory
        error('Cannot find local FCTDdir: %s',FCTDdir);
    end
    if ~exist(dirs.fctd_deployment,'dir') % Bethan 20 June 2021: Added FCTD directory
        error('Cannot find local FCTDdir: %s',FCTDdir);
    end
end

%% rsync remote and local directories
% OR use this option to copy files from one big directory into individual
% deployment directories
if rSync
    % NC - added optional fileStr - find first file to start rsync
    if fileStr
        % List all the files in RawDir and find the first one that matches
        % suffixSearch
        raw_list = dir(fullfile(dirs.raw_incoming,suffixSearch));

        % Identify the file you want
        names = {raw_list.name};
        times = datenum({raw_list.date});
        names(~cellfun('isclass', names, 'char')) = {''};  % Care for non-strings
        matchC = reshape(strfind(names, str_to_match), size(raw_list));
        %Indices of all files that match string:
        match  = ~cellfun('isempty', matchC);
        % Which of these is the earliest file? I think it will always be
        % sorted such that find(match,1,'first') is always the first. If
        % this isn't always the case, use 'times' to sort the matches by
        % date created and find the earliest
        ind = find(match,1,'first');
        first_file = fullfile(dirs.raw_incoming,names{ind});
        %end

        % Rsync the files to a directory specific to the deployment
        if ~exist(dirs.raw_copy,'dir')
            mkdir(dirs.raw_copy,'raw');
        end
        com = sprintf('/usr/bin/rsync -auq --progress --files-from=<(find %s -newer %s -type f -exec basename {} %s) %s %s',dirs.raw_incoming,first_file,'\;',dirs.raw_incoming,dirs.raw_copy);
        fprintf(1,'Running: %s\n',com);
        unix(com);
        fprintf(1,'Done\n');

    end
end


%% Loop through files in the deployment raw directory and convert to mat and fctd_mat
try
    myASCIIfiles = dir(fullfile(dirs.raw_copy, suffixSearch));
    raw_dir = dirs.raw_copy;
catch
    myASCIIfiles = dir(fullfile(dirs.raw_incoming, suffixSearch));
    raw_dir = dirs.raw_incoming;
end

if ~isempty(myASCIIfiles)
    for i=1:length(myASCIIfiles)
        % Sometimes you end up with hidden files if you're copying/pasting files.
        % Don't include those in myASCIIfiles
        if strcmp(myASCIIfiles(i).name(1),'.')
            %do nothing if filename begins with '.'
        else
            indSuffix = strfind(myASCIIfiles(i).name,suffixStr);

            base = myASCIIfiles(i).name(1:indSuffix-1);
            myMATfile = dir(fullfile(dirs.mat, [base '.mat']));

            % If the mat file exists, load it to get the file size of the raw file
            % that created it.
            if ~isempty(myMATfile)
                load(fullfile(myMATfile.folder,myMATfile.name),'raw_file_info')
            end

            % Convert new data to .mat if:
            %   (1) there is no .mat file to match the current raw file, or
            %   (2) the current raw file size is larger than what is saved in mat data 'raw_file_size'
            if ~isempty(myMATfile) && myASCIIfiles(i).bytes<=raw_file_info.bytes && i~=length(myASCIIfiles)
                % If the MAT file exists and its saved raw file is equal to the ascii file, skip recoverting. All
                % of the raw data has already been converted
                % (DO NOTHING.)


            elseif (~isempty(myMATfile) && myASCIIfiles(i).bytes>raw_file_info.bytes) || isempty(myMATfile)
                % If the MAT file exists already but is older than the raw data
                % file, it will be reconverted. OR, if the MAT file does not exist
                % yet, it will be converted.

                % % For debugging
                % debug.base_name{i} = base;
                % debug.rawraw_date(i) = myRAWRAWfile.datenum;
                % debug.rawraw_bytes(i) = myRAWRAWfile.bytes;
                % debug.raw_date(i) = myASCIIfiles(i).datenum;
                % debug.raw_bytes(i) = myASCIIfiles(i).bytes;
                % if ~isempty(myMATfile)
                % debug.mat_date(i) = myMATfile.datenum;
                % debug.mat_bytes(i) = myMATfile.bytes;
                % end
                % debug.conversion_happens(i) = 1;
                % % End for debugging

                %fprintf(1,'Converting %s%s\n',dirs.mat,myMATfile.name);
                fprintf(1,'%s%s --> %s.mat\n', base,suffixStr,base')

                % Read file and save data in matData structure
                filename = fullfile(raw_dir,myASCIIfiles(i).name);
                matData = read_data_file(filename,Meta_Data,version);

                % Add raw_file_info to matData - NC added 8 Aug. 2022
                matData.raw_file_info.bytes = myASCIIfiles(i).bytes;

                % Option to display file data
                if display_file_data_flag
                    display_file_data(myASCIIfiles,i,matData)
                end

                % Save data in .mat file
                save(fullfile(MatDir,base),'-struct','matData')
                fprintf(1,'%s: Wrote  %s.mat\n',datestr(now,'YY.mm.dd HH:MM:SS'),base);

                %Empty contents of matData structure
                use matData

                % Calculate microstructure data - NC added 16 May 2022
                if calc_micro
                    matData.micro = epsiProcess_calc_turbulence(Meta_Data,matData,0);
                    save(fullfile(MatDir,base),'-struct','matData')
                end

                % Update the .mat file time index
                if ~isempty(epsi) && isfield(epsi,'time_s')
                    epsiProcess_update_TimeIndex(MatDir,base,epsi);
                elseif isempty(epsi) &&  ~isempty(ctd) && isfield(ctd,'time_s') %For the case where the is no epsi data, but there is ctd data
                    epsiProcess_update_TimeIndex(MatDir,base,ctd);
                end

                % Update pressure timeseries
                if ~isempty(ctd) && isfield(ctd,'dnum')
                    epsiProcess_update_PressureTimeseries(Meta_Data,MatDir,ctd,Meta_Data.paths.profiles)
                end

                % Save file for FCTD Format %%%%%% (Bethan 20 June 2021)
                if doFCTD && ~isempty(ctd) && isfield(ctd,'time_s')
                    % Make FCTD .mat data and save it in cruise-long fctd
                    % directory
                    FCTD = make_FCTD_mat(matData,dirs.fctd_cruise,base,cruise_specifics);

                    % Also save it in deployment fctd directory
                    save(fullfile(dirs.fctd_deployment,base),'FCTD');
                end %end if doFCTD

            end %end if the data should be converted
        end %end if name of file does not begin with '.'
    end %end loop through files

    %eval(['save ' fullfile('~/Desktop',strrep(strrep(datestr(now),':','_'),' ','_')) ' debug'])
end %end if there are files

end



% ------------------------------
% ----- SUBFUNCTIONS -----------
% ------------------------------
function  [matData] = read_data_file(filename,Meta_Data,version)

switch version
    case 4
        [matData] = mod_som_read_epsi_files_v4(filename,Meta_Data);
        use matData
    case 3
        [epsi,ctd,alt,act,vnav,gps] = mod_som_read_epsi_files_v3(filename,Meta_Data);
        matData.epsi = epsi;
        matData.ctd  = ctd;
        matData.alt  = alt;
        matData.act  = act;
        matData.vnav = vnav;
        matData.gps  = gps;
    case 2
        [epsi,ctd,alt,act]=mod_som_read_epsi_files_v2(filename,Meta_Data);
        matData.epsi = epsi;
        matData.ctd  = ctd;
        matData.alt  = alt;
        matData.act  = act;
        vnav = [];
        gps = [];
    case 1
        [epsi,ctd,alt,act]=mod_som_read_epsi_files_v1(filename,Meta_Data);
        matData.epsi = epsi;
        matData.ctd = ctd;
        matData.alt = alt;
        matData.act = act;
        vnav = [];
        gps = [];
    case 0
        EPSI = mod_read_epsi_raw(filename,Meta_Data);
        % Save as standalone structures too, for saving individual
        % files in following steps
        epsi = EPSI.epsi;
        ctd = EPSI.aux1;
        act = [];
        alt = [];
        vnav = [];
        gps = [];

        t0 = Meta_Data.starttime;
        epsi.time_s=epsi.EPSInbsample/Meta_Data.PROCESS.Fs_epsi;
        ctd.time_s=ctd.Aux1Stamp/Meta_Data.PROCESS.Fs_epsi;

        if (t0==0)
            epsi.dnum = epsi.time;
            ctd.dnum = ctd.time;
        else
            ctd.dnum = ctd.time_s/86400 +t0;
            epsi.dnum = epsi.time_s/86400 +t0;
        end

        % Add extra ctd variables
        % Define a constant for salinity calculation
        c3515 = 42.914;
        ctd.S    = sw_salt(ctd.C*10./c3515,ctd.T,ctd.P);
        ctd.th   = sw_ptmp(ctd.S,ctd.T,ctd.P,0);
        ctd.sgth  = sw_pden(ctd.S,ctd.T,ctd.P,0);
        ctd.dPdt = [0; diff(ctd.P)./diff(ctd.time_s)];

        % NC 17 July 2021 - added ctd.z and ctd.dzdt.
        % get_scan_spectra.m will use dzdt to define fall speed w.
        if ~isfield(Meta_Data.PROCESS,'latitude')
            warning('Need latitude to get depth from pressure data. Add to MetaProcess text file. Pursue without')
            ctd.z    = ctd.P;
            ctd.dzdt = [0; diff(ctd.z)./diff(ctd.time_s)];
        else
            ctd.z    = sw_dpth(ctd.P,Meta_Data.PROCESS.latitude);
            ctd.dzdt = [0; diff(ctd.z)./diff(ctd.time_s)];
        end

        matData.epsi = epsi;
        matData.ctd = ctd;
end
end %end read_data_file
% ---------------------------------



function  [] = display_file_data(myASCIIfiles,i,matData)

use matData % Empty contents of matData structure

% Display file size, time, pressure, altimeter
try
    disp('~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ')
    disp(['FILE SIZE: ' num2str(myASCIIfiles(i).bytes)])
    disp(['TIME: ' datestr(ctd.dnum(end))])
    disp(['PRESSURE: ' num2str(ctd.P(end))])
catch
end
try
    disp(['ALTIMETER: ' num2str(alt.dst(end))])
catch
end
end %end display_file_data
% ------------------------------------


% function [FCTD] =  make_FCTD_mat(matData,FCTDdir,base,cruise_specifics)
% %THIS IS NOW A SEPARATE FUNCTION
%
% use matData %Empty contents of matData structure
%
% % Get CTD data
% FCTD.time=ctd.dnum;
% FCTD.pressure=ctd.P;
% FCTD.temperature=ctd.T;
% FCTD.conductivity=ctd.C;
%
% % Get altimeter data
% if ~isempty(alt) && isfield(alt,'time_s')
%     FCTD.altDist=interp1(alt.dnum,alt.dst,ctd.dnum);
% else
%     FCTD.altTime=nan(length(ctd.dnum),1);
% end
%
% % Add VectorNav data
% if ~isempty(vnav) && isfield(vnav,'time_s')
%     diff_not_neg = [0;diff(vnav.dnum)]>0;
%     keep = ~isnan(vnav.dnum) & ~isinf(vnav.dnum) & diff_not_neg;
%     for ix=1:3
%         FCTD.compass(:,ix)=interp1(vnav.dnum(keep),vnav.compass(keep,ix),ctd.dnum);
%         FCTD.gyro(:,ix)=interp1(vnav.dnum(keep),vnav.gyro(keep,ix),ctd.dnum);
%         FCTD.acceleration(:,ix)=interp1(vnav.dnum(keep),vnav.acceleration(keep,ix),ctd.dnum)./9.81;
%     end
% else
%     FCTD.gyro=nan(length(ctd.dnum),3);
%     FCTD.acceleration=nan(length(ctd.dnum),3);
%     FCTD.compass=nan(length(ctd.dnum),3);
% end
%
% % Add GPS data
% if ~isempty(gps) && isfield(gps,'gpstime')
%     FCTD.GPS.longitude=interp1(gps.gpstime,gps.longitude,ctd.dnum);
%     FCTD.GPS.latitude=interp1(gps.gpstime,gps.latitude,ctd.dnum);
% else
%     FCTD.GPS.longitude=nan(length(ctd.dnum),1);
%     FCTD.GPS.latitude=nan(length(ctd.dnum),1);
% end
%
% % Extra outputs for BLT 2021 cruise
% if cruise_specifics.blt_2021
%     % Microconductivity and Fluorometer
%     %
%     % On BLT 2021, microconductivity sensor was on shear
%     % channel 2 of epsi and fluorometer was on shear channel 1. This step interpolates that data to
%     % the same time array as the rest of the data, but since it
%     % has a 20x faster sampling rate than the SBE (320 Hz vs 16
%     % Hz), it actually becomes and N x 20 array - there are 20
%     % uConductivity/fluorometer data points for every 1 SBE data point. We
%     % also save time_fast as an N x 20 array.
%     time_fast = linspace(ctd.dnum(1),ctd.dnum(end),length(ctd.dnum)*20);
%     FCTD.time_fast = time_fast;
%
%     % Interpolate data that is not nan, not inf, and where time
%     % is increasing
%     diff_not_neg = [0;diff(epsi.dnum)]>0;
%     keep = ~isnan(epsi.dnum) & ~isinf(epsi.dnum) & diff_not_neg;
%
%     if ~isempty(epsi) && isfield(epsi,'s2_count') && ~isempty(ctd)
%         FCTD.uConductivity=reshape(interp1(epsi.dnum(keep),double(epsi.s2_count(keep)),time_fast),20,[])';
%     else
%         FCTD.uConductivity=nan(length(ctd.dnum),20);
%         disp(['No uConductivity data ' myASCIIfiles(i).name]);
%     end
%
%     if ~isempty(epsi) && isfield(epsi,'s1_volt')  && ~isempty(ctd)
%         FCTD.fluorometer=reshape(interp1(epsi.dnum(keep),epsi.s1_volt(keep),time_fast),20,[])';
%     else
%         FCTD.fluorometer=nan(length(ctd.dnum),20);
%         disp(['No fluorometer data ' myASCIIfiles(i).name]);
%     end
% end %end if blt_2021
%
% % Save FCTD mat files to the new FCTD mat directory FCTDmat
% myFCTDMATfile = fullfile(FCTDdir,base);
% save(myFCTDMATfile,'FCTD');
% fprintf(1,'%s: Wrote  %s%s\n\n',datestr(now,'YY.mm.dd HH:MM:SS'), FCTDdir,myFCTDMATfile);
%
% % Update FCTD .mat time index
% FastCTD_UpdateMATFileTimeIndex(FCTDdir,base,FCTD);
%
% end %end make_FCTD_mat
% % ---------------------------------



% function FastCTD_UpdateMATFileTimeIndex(dirname,filename,FCTD)
% if exist([dirname '/FastCTD_MATfile_TimeIndex.mat'],'file')
%     load([dirname '/FastCTD_MATfile_TimeIndex.mat']);
%     ind = strncmp(filename,FastCTD_MATfile_TimeIndex.filenames,length(filename));
%     if sum(ind) ~= 1
%         FastCTD_MATfile_TimeIndex.filenames = [FastCTD_MATfile_TimeIndex.filenames; {filename}];
%         FastCTD_MATfile_TimeIndex.timeStart = cat(1,FastCTD_MATfile_TimeIndex.timeStart,FCTD.time(1));
%         FastCTD_MATfile_TimeIndex.timeEnd = cat(1,FastCTD_MATfile_TimeIndex.timeEnd,FCTD.time(end));
%     else
%         FastCTD_MATfile_TimeIndex.timeStart(ind) = FCTD.time(1);
%         FastCTD_MATfile_TimeIndex.timeEnd(ind) = FCTD.time(end);
%     end
% else
%     FastCTD_MATfile_TimeIndex.filenames = {filename};
%     FastCTD_MATfile_TimeIndex.timeStart = FCTD.time(1);
%     FastCTD_MATfile_TimeIndex.timeEnd = FCTD.time(end);
% end
% save([dirname '/FastCTD_MATfile_TimeIndex.mat'],'FastCTD_MATfile_TimeIndex');
% end %end FastCTD_UpdateMATFileTimeIndex
% % ------------------------------------
