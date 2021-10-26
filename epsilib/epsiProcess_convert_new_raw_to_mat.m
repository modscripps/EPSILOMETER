function [matData] = epsiProcess_convert_new_raw_to_mat(dirs,Meta_Data,varargin)
% epsiProcess_convert_new_raw_to_mat
%
% Nicole Couto adapted from FCTD_MakeMatFromRaw.m
% May-July 2021
% --------------------------------------------------------
% INPUTS:
%   dirs = {RawDir;
%           RawDirAway;
%           MatDir;
%           FctdMatDir;}
%   Meta_Data: epsi Meta_Data structure
%
% OPTIONAL INPUTS:
%   'noSync' - don't rsync files from RawDir to RawDir Away
%   'noGrid' - don't grid fctd data (To do: remove this. We don't use this
%              anymore)
%   'doFCTD' - make FCTD files compatible with MOD FCTD processing library
%   'fileStr', 'string_to_match' - use this flag along with a search string
%       if you want to limit rsyncing files to only those that match a certain
%       pattern. This is useful if you're storing all the raw data from a
%       cruise in RawDir, but you want to separate deployments into different
%       RawDirAway. You could specify only files from July 17th as
%       'fileStr','*07_17*
%
% Written by Jody Klymak
% Updated 2011 06 21 by San Nguyen
% Updated 2012 09 29 by San Nguyen for EquatorMix2012

% NC - Make matData for output even if there is no new data
matData.epsi = [];
matData.ctd = [];
matData.alt = [];
matData.act = [];
matData.vnav = [];
matData.gps = [];

% NC - Only rsync files with the desired suffix
suffixStr = Meta_Data.rawfileSuffix; %ex. *.raw, *.ascii, etc
suffixSearch = ['*' suffixStr];

% NC - make sure all dirs have a / at the end
for ii=1:numel(dirs)
    dirs{ii} = strrep([dirs{ii},'/'],'//','/');
end

persistent argsNameToCheck;
if isempty(argsNameToCheck);
    argsNameToCheck = {'noSync','noGrid','doFCTD','fileStr','version'}; %Add flag for doFCTD - makes mat files in FCTD processing format
end

% TODO PLease comment on these parameters
rSync = true;
doGrid = false;
doFCTD = false;
version = 3;

index = 1;
n_items = nargin-2;

% Loop through the number of varargin arguments, check which
% argsNameToCheck it matches, and switch on/off rSync,doGrid, and doFCTD
% accordingly
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
        case 2 % noGrid
            doGrid = false;
            index = index +1;
            n_items = n_items-1;
        case 3 %doFCTD
            doFCTD = true;
            index = index+1;
            n_items = n_items-1;
        case 4 %fileStr % NC - added optional fileStr input
            fileStr = true;
            % Find the index of varargin that = 'fileStr'. The following
            % index contains 'str_to_match'
            idxFlag = find(cell2mat(cellfun(@(C) ~isempty(strfind(C,'fileStr')),varargin,'uniformoutput',0)));
            str_to_match = varargin{idxFlag+1};
            index = index+2; %+2 because the following varargin will str_to_match
            n_items = n_items-2;
        case 5 %version % NC 10/11/21 - added optional version input
            % Find the index of varargin that = 'version'. The following
            % index contains the version number
            idxFlag = find(cell2mat(cellfun(@(C) ~isempty(strfind(C,'version')),varargin,'uniformoutput',0)));
            version = varargin{idxFlag+1};
            index = index+2; %+2 because the following varargin will be the version number
            n_items = n_items-2;
    end
end

if rSync
    RawDir = [dirs{1}];
    RawDirDuplicate = [dirs{2}];
    MatDir = [dirs{3}];
    if doFCTD
        FCTDdir = [dirs{4}]; % Bethan 20 June 2021: Added FCTD directory
    end
    %     if doGrid %NC - doGrid is an option in FCTD processing, but we don't
    %     need it at this step of Epsi processing
    %        GridDir = [dirs{5}]; %Bethan 20 June 2021: Changed from dirs{4} to dirs{5} to account for FCTD
    %     end
    
    % check for valid directories
    if ~exist(RawDirDuplicate,'dir')
        error('Cannot find remote dir: %s',RawDirDuplicate);
    end
else
    RawDir = [dirs{1}];
    MatDir = [dirs{2}];
    if doFCTD
        FCTDdir=[dirs{3}]; % Bethan 20 June 2021: Added FCTD directory
    end
    %     if doGrid %NC - doGrid is an option in FCTD processing, but we don't
    %     need it at this step of Epsi processing
    %         GridDir = [dirs{4}]; %Bethan 20 June 2021: Changed from dirs{3} to dirs{4} to account for FCTD
    %     end
end

if ~exist(RawDir,'dir')
    error('Cannot find local RawDir: %s',RawDir);
end
if ~exist(MatDir,'dir')
    error('Cannot find local MatDir: %s',MatDir);
end
if doGrid && ~exist(GridDir,'dir')
    error('Cannot find local GridDir: %s',GridDir);
end
if doFCTD && ~exist(FCTDdir,'dir')
    error('Cannot find local FCTDdir: %s',FCTDdir);
end

% rsync remote and local directories.  You must connect to the server
% to make this work.
if rSync
    % NC - added optional fileStr input
    if fileStr
        suffixSearch = strrep([str_to_match,suffixSearch],'**','*'); %If you ended up with side-by-side *, delete one
    end
    com = sprintf('/usr/bin/rsync -av  --include ''%s'' --exclude ''*'' %s %s',suffixSearch,RawDir,RawDirDuplicate);  %NC - rsync only the files ending in .raw
    
    fprintf(1,'Running: %s\n',com);
    unix(com);
    fprintf(1,'Done\n');
    RawDir = RawDirDuplicate;
end

myASCIIfiles = dir([RawDir, suffixSearch]);

for i=1:length(myASCIIfiles)
    indSuffix = strfind(myASCIIfiles(i).name,suffixStr);
    
    base = myASCIIfiles(i).name(1:indSuffix-1);
    myMATfile = dir([MatDir base '.mat']);
    if doFCTD
        myFCTDMATfile = dir([FCTDdir base '.mat']);
    end
    
    % if the MAT files are older than the data files, they will be retranslated
    if (~isempty(myMATfile) && datenum(myASCIIfiles(i).date)>datenum(myMATfile.date))
        fprintf(1,'Retranslating %s%s\n',MatDir,myMATfile.name);
        
        filename = fullfile(RawDir,myASCIIfiles(i).name);
        [matData, epsi,ctd,alt,act,vnav,gps] = read_data_file(filename,Meta_Data,version);
        
        %             % Display file size, time, pressure, altimeter
        %             try
        %                 disp('~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ')
        %                 disp(['FILE SIZE: ' num2str(myASCIIfiles(i).bytes)])
        %                 disp(['TIME: ' datestr(ctd.dnum(end))])
        %                 disp(['PRESSURE: ' num2str(ctd.P(end))])
        %             catch
        %             end
        %             try
        %                 disp(['ALTIMETER: ' num2str(alt.dst(end))])
        %             catch
        %             end
        
        if ~isempty(epsi) && isfield(epsi,'time_s')
            save([MatDir  base '.mat'],'epsi','ctd','alt');
            epsiProcess_update_TimeIndex(MatDir,base,epsi);
            fprintf(1,'%s: Wrote  %s%s.mat\n\n',datestr(now,'YY.mm.dd HH:MM:SS'),MatDir,base);
        elseif isempty(epsi) &&  ~isempty(ctd) && isfield(ctd,'time_s') %For the case where the is no epsi data, but there is ctd data
            save([MatDir  base '.mat'],'epsi','ctd','alt','vnav');
            epsiProcess_update_TimeIndex(MatDir,base,ctd);
        end
        
        % Update pressure timeseries
        if ~isempty(ctd) && isfield(ctd,'dnum')
            epsiProcess_update_PressureTimeseries(MatDir,ctd)
        end
        
        %%%%% Save files for FCTD Format %%%%%% (Bethan 20 June 2021)
        if doFCTD
            time_offset = 0; % 2021 07 03 SAN added to correct for time for the current deployment
            if ~isempty(ctd) && isfield(ctd,'time_s')
                %Create new structure FCTD with data from alt, ctd and epsi renamed to work with current FCTD processing
                
                % Get CTD data
                FCTD.time=ctd.dnum+time_offset; %Currently just ctdtime (seconds since powered on - will need to change this once we have gps data)
                FCTD.pressure=ctd.P;
                FCTD.temperature=ctd.T;
                FCTD.conductivity=ctd.C;
                
                % Get altimeter data
                if ~isempty(alt) && isfield(alt,'alttime')
                    FCTD.altDist=interp1(alt.altdnum,alt.dst,ctd.dnum);
                else
                    FCTD.altTime=nan(length(ctd.dnum),1);
                    disp(['No alt data ' myASCIIfiles(i).name]);
                end
                
                % Get microconductivity (this is saved in shear channel 2
                % of epsi - needs to be interpolated onto the same time
                % base (x20 to account for higher sampling rate) as the
                % rest of the data
                % THIS IS SAVED IN VOLTS NOT IN MICROCONDUCTIVITY UNITS SO
                % WILL NEED TO MAKE SURE PROCESSING FURTHER DOWN THE LINE
                % ACCOUNTS FOR THIS
                if ~isempty(epsi) && isfield(epsi,'s2_count') && ~isempty(ctd)
                    ucontime=linspace(ctd.dnum(1),ctd.dnum(end),length(ctd.dnum)*20);
                    FCTD.uConductivity=reshape(interp1(epsi.dnum,double(epsi.s2_count),ucontime),20,[])';
                    clear ucontime
                else
                    FCTD.uConductivity=nan(length(ctd.dnum),20);
                    disp(['No uConductivity data ' myASCIIfiles(i).name]);
                end
                
                
                % If we want the fluorometer to be outputted as well then
                % that is saved in shear channel 1
                % Currently saving in the same format as uCond but we may
                % not need it to be so high resolution
                if ~isempty(epsi) && isfield(epsi,'s1_volt')  && ~isempty(ctd)
                    fluortime=linspace(ctd.dnum(1),ctd.dnum(end),length(ctd.dnum)*20);
                    FCTD.fluorometer=reshape(interp1(epsi.dnum,epsi.s1_volt,fluortime),20,[])';
                    clear fluortime
                else
                    FCTD.fluorometer=nan(length(ctd.dnum),20);
                    disp(['No fluorometer data ' myASCIIfiles(i).name]);
                end
                
                % Get VectorNav data
                
                % Add vnav.vnavdnum data, interpolate to CTD
                if ~isempty(vnav) && isfield(vnav,'vnavtime')
                    for ix=1:3
                        FCTD.compass(:,ix)=interp1(vnav.vnavdnum,vnav.compass(:,ix),ctd.dnum);
                        FCTD.gyro(:,ix)=interp1(vnav.vnavdnum,vnav.gyro(:,ix),ctd.dnum);
                        FCTD.acceleration(:,ix)=(interp1(vnav.vnavdnum,vnav.acceleration(:,ix),ctd.dnum))./9.81;
                    end
                else
                    FCTD.gyro=nan(length(ctd.dnum),3);
                    FCTD.acceleration=nan(length(ctd.dnum),3);
                    FCTD.compass=nan(length(ctd.dnum),3);
                end
                
                %%%%%% ADD GPS WHEN WE HAVE THAT DATA %%%%%%%
                
                if ~isempty(vnav) && isfield(vnav,'vnavtime')
                    for ix=1:3
                        FCTD.compass(:,ix)=interp1(vnav.vnavdnum,vnav.compass(:,ix),ctd.dnum);
                        FCTD.gyro(:,ix)=interp1(vnav.vnavdnum,vnav.gyro(:,ix),ctd.dnum);
                        FCTD.acceleration(:,ix)=(interp1(vnav.vnavdnum,vnav.acceleration(:,ix),ctd.dnum))./9.81;
                    end
                else
                    FCTD.gyro=nan(length(ctd.dnum),3);
                    FCTD.acceleration=nan(length(ctd.dnum),3);
                    FCTD.compass=nan(length(ctd.dnum),3);
                end
                
                % Add GPS data
                
                if ~isempty(gps) && isfield(gps,'gpstime')
                    FCTD.GPS.longitude=interp1(gps.gpstime,gps.longitude,ctd.dnum);
                    FCTD.GPS.latitude=interp1(gps.gpstime,gps.latitude,ctd.dnum);
                else
                    FCTD.GPS.longitude=nan(length(ctd.dnum),1);
                    FCTD.GPS.latitude=nan(length(ctd.dnum),1);
                end
                
                
                % Save FCTD mat files to the new FCTD mat directory FCTDmat
                save([FCTDdir  base '.mat'],'FCTD');
                FastCTD_UpdateMATFileTimeIndex(FCTDdir,base,FCTD);
                fprintf(1,'%s: Wrote  %s%s\n\n',datestr(now,'YY.mm.dd HH:MM:SS'), FCTDdir,myFCTDMATfile.name);
                %                 if doGrid
                %                     FCTD_GridData = Epsi_GridData(FCTD);
                %                     save([GridDir base '.mat'],'FCTD_GridData');
                %                     epsiProcess_update_TimeIndex(GridDir,base,FCTD_GridData);
                %                     fprintf(1,'%s: Wrote  %s%s.mat\n\n',datestr(now,'YY.mm.dd HH:MM:SS'),GridDir,base);
                %                 end
            end %end Bethan's addition for FCTD data
        end %end if doFCTD
        
        
        
        %         catch err
        %             disp(['So... this is the error for retranlating file ' myASCIIfiles(i).name]);
        %             disp(err);
        %             for j = 1:length(err.stack)
        %                 disp([num2str(j) ' ' err.stack(j).name ' ' num2str(err.stack(j).line)]);
        %             end
        %             error('There was an error. See stack above')
        %         end
        
        % If the files are new then a new MAT file will be created
    elseif isempty(myMATfile)
        fprintf(1,'Translating %s%s\n',RawDir,myASCIIfiles(i).name);
        
        filename = fullfile(RawDir,myASCIIfiles(i).name);
        [matData, epsi,ctd,alt,act,vnav,gps] = read_data_file(filename,Meta_Data,version);
        

        
        %             % Display file size, time, pressure, altimeter
        %             try
        %                 disp('~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ')
        %                 disp(['FILE SIZE: ' num2str(myASCIIfiles(i).bytes)])
        %                 disp(['TIME: ' datestr(ctd.dnum(end))])
        %                 disp(['PRESSURE: ' num2str(ctd.P(end))])
        %             catch
        %             end
        %             try
        %                 disp(['ALTIMETER: ' num2str(alt.dst(end))])
        %             catch
        %             end
        
        % Save the data in a mat file with the same name as the raw file
        if ~isempty(epsi) && isfield(epsi,'time_s')
            save([MatDir  base '.mat'],'epsi','ctd','alt','vnav','gps');
            epsiProcess_update_TimeIndex(MatDir,base,epsi);
            fprintf(1,'%s: Wrote  %s%s.mat\n\n',datestr(now,'YY.mm.dd HH:MM:SS'),MatDir,base);
        elseif isempty(epsi) &&  ~isempty(ctd) && isfield(ctd,'time_s') %For the case where the is no epsi data, but there is ctd data
            save([MatDir  base '.mat'],'epsi','ctd','alt','vnav','gps');
            epsiProcess_update_TimeIndex(MatDir,base,ctd);
            fprintf(1,'%s: Wrote  %s%s.mat\n\n',datestr(now,'YY.mm.dd HH:MM:SS'),MatDir,base);
        end
        
        % Update pressure timeseries
        if ~isempty(ctd) && isfield(ctd,'dnum')
            epsiProcess_update_PressureTimeseries(MatDir,ctd)
        end
        
        
        if doFCTD
            % Bethan's addition for FCTD data
            time_offset = 0; % 2021 07 03 SAN added to correct for time for the current deployment
            if ~isempty(ctd) && isfield(ctd,'time_s')
                %Create new structure FCTD with data from alt, ctd and epsi renamed to work with current FCTD processing
                
                % Get CTD data
                FCTD.time=ctd.dnum+time_offset;
                FCTD.pressure=ctd.P;
                FCTD.temperature=ctd.T;
                FCTD.conductivity=ctd.C;
                
                % Get altimeter data
                if ~isempty(alt) && isfield(alt,'time_s')
                    FCTD.altDist=interp1(alt.dnum,alt.dst,ctd.dnum);
                else
                    FCTD.altTime=nan(length(ctd.dnum),1);
                    disp(['No alt data ' myASCIIfiles(i).name]);
                end
                
                % Get microconductivity (this is saved in shear channel 2
                % of epsi - needs to be interpolated onto the same time
                % base (x20 to account for higher sampling rate) as the
                % rest of the data
                % THIS IS SAVED IN VOLTS NOT IN MICROCONDUCTIVITY UNITS SO
                % WILL NEED TO MAKE SURE PROCESSING FURTHER DOWN THE LINE
                % ACCOUNTS FOR THIS
                if ~isempty(epsi) && isfield(epsi,'s2_count')  && ~isempty(ctd)
                    ucontime=linspace(ctd.dnum(1),ctd.dnum(end),length(ctd.dnum)*20);
                    FCTD.uConductivity=reshape(interp1(epsi.dnum,double(epsi.s2_count),ucontime),20,[])';
                    clear ucontime
                else
                    FCTD.uConductivity=nan(length(ctd.dnum),20);
                    disp(['No uConductivity data ' myASCIIfiles(i).name]);
                end
                
                
                % If we want the fluorometer to be outputted as well then
                % that is saved in shear channel 1
                % Currently saving in the same format as uCond but we may
                % not need it to be so high resolution
                if ~isempty(epsi) && isfield(epsi,'s1_volt')  && ~isempty(ctd)
                    fluortime=linspace(ctd.dnum(1),ctd.dnum(end),length(ctd.dnum)*20);
                    FCTD.fluorometer=reshape(interp1(epsi.dnum,epsi.s1_volt,fluortime),20,[])';
                    clear fluortime
                else
                    FCTD.fluorometer=nan(length(ctd.dnum),20);
                    disp(['No fluorometer data ' myASCIIfiles(i).name]);
                end
                
                % Get VectorNav data
                
                
                %Bethan to do:
                % Add vnav.vnavdnum data, interpolate to CTD
                if ~isempty(vnav) && isfield(vnav,'time_s')
                    for ix=1:3
                        FCTD.compass(:,ix)=interp1(vnav.dnum,vnav.compass(:,ix),ctd.dnum);
                        FCTD.gyro(:,ix)=interp1(vnav.dnum,vnav.gyro(:,ix),ctd.dnum);
                        FCTD.acceleration(:,ix)=(interp1(vnav.dnum,vnav.acceleration(:,ix),ctd.dnum))./9.81;
                    end
                else
                    FCTD.gyro=nan(length(ctd.dnum),3);
                    FCTD.acceleration=nan(length(ctd.dnum),3);
                    FCTD.compass=nan(length(ctd.dnum),3);
                end
                
                % Add GPS data
                
                if ~isempty(gps) && isfield(gps,'dnum')
                    FCTD.GPS.longitude=interp1(gps.dnum,gps.longitude,ctd.dnum);
                    FCTD.GPS.latitude=interp1(gps.dnum,gps.latitude,ctd.dnum);
                else
                    FCTD.GPS.longitude=nan(length(ctd.dnum),1);
                    FCTD.GPS.latitude=nan(length(ctd.dnum),1);
                end
                
                
                % Save FCTD mat files to the new FCTD mat directory FCTDmat
                save([FCTDdir  base '.mat'],'FCTD');
                FastCTD_UpdateMATFileTimeIndex(FCTDdir,base,FCTD);
                fprintf(1,'%s: Wrote  %s%s\n\n',datestr(now,'YY.mm.dd HH:MM:SS'), FCTDdir,myFCTDMATfile.name);
                %                 if doGrid
                %                     FCTD_GridData = Epsi_GridData(FCTD);
                %                     save([GridDir base '.mat'],'FCTD_GridData');
                %                     epsiProcess_update_TimeIndex(GridDir,base,FCTD_GridData);
                %                     fprintf(1,'%s: Wrote  %s%s.mat\n\n',datestr(now,'YY.mm.dd HH:MM:SS'),GridDir,base);
                %                 end
                
            end %end Bethan's addition for FCTD data
        end %end if doFCTD
        
        %         catch err
        %             disp(['So... this is the error for tranlating file ' myASCIIfiles(i).name]);
        %             disp(err);
        %             for j = 1:length(err.stack)
        %                 disp([num2str(j) ' ' err.stack(j).name ' ' num2str(err.stack(j).line)]);
        %             end
        %             error('There was an error. See stack above.')
        %         end
    end
    clear FCTD
end

clear epsi ctd alt act vnav;
end



% ------------------------------
% ----- SUBFUNCTIONS -----------
% ------------------------------
function  [matData, epsi,ctd,alt,act,vnav,gps] = read_data_file(filename,Meta_Data,version)

switch version
    case 3
        [epsi,ctd,alt,act,vnav,gps] = mod_som_read_epsi_files_v3(filename,Meta_Data);
        matData.epsi = epsi;
        matData.ctd = ctd;
        matData.alt = alt;
        matData.act = act;
        matData.vnav = vnav;
        matData.gps = gps;
    case 2
        [epsi,ctd,alt,act]=mod_som_read_epsi_files_v2(filename,Meta_Data);
        matData.epsi = epsi;
        matData.ctd = ctd;
        matData.alt = alt;
        matData.act = act;
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
        
        
        epsi.dnum = epsi.time;
        t0 = Meta_Data.starttime;
        epsi.time_s = (epsi.dnum-t0)*(24*60*60);
        
        ctd.dnum = ctd.time;
        ctd.time_s = (ctd.dnum-t0)*(24*60*60);
        
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
            error('Need latitude to get depth from pressure data. Add to MetaProcess text file.')
        else
            ctd.z    = sw_dpth(ctd.P,Meta_Data.PROCESS.latitude);
            ctd.dzdt = [0; diff(ctd.z)./diff(ctd.time_s)];
        end
        
        matData.epsi = epsi;
        matData.ctd = ctd;
end
end %end read_data_file
% ---------------------------------



function FastCTD_UpdateMATFileTimeIndex(dirname,filename,FCTD)
if exist([dirname '/FastCTD_MATfile_TimeIndex.mat'],'file')
    load([dirname '/FastCTD_MATfile_TimeIndex.mat']);
    ind = strncmp(filename,FastCTD_MATfile_TimeIndex.filenames,length(filename));
    if sum(ind) ~= 1
        FastCTD_MATfile_TimeIndex.filenames = [FastCTD_MATfile_TimeIndex.filenames; {filename}];
        FastCTD_MATfile_TimeIndex.timeStart = cat(1,FastCTD_MATfile_TimeIndex.timeStart,FCTD.time(1));
        FastCTD_MATfile_TimeIndex.timeEnd = cat(1,FastCTD_MATfile_TimeIndex.timeEnd,FCTD.time(end));
    else
        FastCTD_MATfile_TimeIndex.timeStart(ind) = FCTD.time(1);
        FastCTD_MATfile_TimeIndex.timeEnd(ind) = FCTD.time(end);
    end
else
    FastCTD_MATfile_TimeIndex.filenames = {filename};
    FastCTD_MATfile_TimeIndex.timeStart = FCTD.time(1);
    FastCTD_MATfile_TimeIndex.timeEnd = FCTD.time(end);
end
save([dirname '/FastCTD_MATfile_TimeIndex.mat'],'FastCTD_MATfile_TimeIndex');
end
