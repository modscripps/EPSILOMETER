function Epsi_MakeMatFromRawGPS_andPlot(dirs,Meta_Data,varargin)
% Epsi_MakeMatFromRaw
%
% Nicole Couto adapted from Epsi_MakeMatFromRaw.m
% May 2021
% --------------------------------------------------------
% Epsi_MakeMatFromRaw(DIRS)
% Script to make FCTD "raw" (ASCII) matfiles
% Run: Epsi_MakeMatFromRaw(DIRS)
% where DIRS should be
% DIRS = {RawDir;
%         RawDirAway;
%         MatDir;
%         GridDir;}
%
% Written by Jody Klymak
% Updated 2011 06 21 by San Nguyen
% Updated 2012 09 29 by San Nguyen for EquatorMix2012

% NC - Only rsync files with the desired suffix
suffixStr = Meta_Data.rawfileSuffix; %ex. *.raw, *.ascii, etc
suffixSearch = ['*' suffixStr];

% NC - make sure all dirs have a / at the end
for ii=1:numel(dirs)
    dirs{ii} = strrep([dirs{ii},'/'],'//','/');
end

persistent argsNameToCheck;
if isempty(argsNameToCheck);
    argsNameToCheck = {'noSync','noGrid','doFCTD'}; %Add flag for doFCTD - makes mat files in FCTD processing format
end

rSync = true;
doGrid = false;
doFCTD = false;

index = 1;
n_items = nargin-2;

% Loop through the number of varargin arguments, check which
% argsNameToCheck it matches, and switch on/off rSync,doGrid, and doFCTD
% accordingly
while (n_items > 0)
    argsMatch = strcmpi(varargin{index},argsNameToCheck);
    i = find(argsMatch,1);
    if isempty(i)
        error('MATLAB:Epsi_MakeMatFromRaw:wrongOption','Incorrect option specified: %s', varargin{index});
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
    %com = sprintf('/usr/bin/rsync -av %s %s',RawDir,RawDirDuplicate);
    %com = sprintf('/usr/bin/rsync -av  --include ''%s'' --exclude ''*'' %s %s',suffixSearch,RawDir,RawDirDuplicate);  %The exclude was useful before... but now it actually excludes everything. Leaving this commented in case I need it again
    %com = sprintf('/usr/bin/rsync -av  --include ''%s'' --exclude ''*'' %s %s',suffixSearch,RawDir,RawDirDuplicate);  %The exclude was useful before... but now it actually excludes everything. Leaving this commented in case I need it again
    com = sprintf('/usr/bin/rsync -av  --include ''%s'' %s %s',suffixSearch,RawDir,RawDirDuplicate);  %NC - rsync only the files with the str_to_search and suffix
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
        try
            disp([RawDir myASCIIfiles(i).name]);

            [epsi,ctd,alt,act,vnav,GPS] = mod_som_read_epsi_files_vGPS([RawDir myASCIIfiles(i).name],Meta_Data);
            % Add fileNum

            % Display file size, time, pressure, altimeter
            try
                disp('~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ')
                disp(['FILE SIZE: ' num2str(myASCIIfiles(i).bytes)])
                disp(['TIME: ' datestr(ctd.ctddnum(end))])
                disp(['PRESSURE: ' num2str(ctd.P(end))])
            catch
            end
            try
                disp(['ALTIMETER: ' num2str(alt.dst(end))])
            catch
            end

            % Update plot



            if ~isempty(epsi) && isfield(epsi,'epsitime')
                save([MatDir  base '.mat'],'epsi','ctd','alt');
                Epsi_UpdateMATFileTimeIndex(MatDir,base,epsi.epsitime);
                fprintf(1,'%s: Wrote  %s%s.mat\n\n',datestr(now,'YY.mm.dd HH:MM:SS'),MatDir,base);
            elseif isempty(epsi) &&  ~isempty(ctd) && isfield(ctd,'ctdtime') %For the case where the is no epsi data, but there is ctd data
                save([MatDir  base '.mat'],'epsi','ctd','alt','vnav');
                Epsi_UpdateMATFileTimeIndex(MatDir,base,ctd.ctdtime);
            end


            %%%%% Save files for FCTD Format %%%%%% (Bethan 20 June 2021)
            if doFCTD
                if ~isempty(ctd) && isfield(ctd,'ctdtime')
                    %Create new structure FCTD with data from alt, ctd and epsi renamed to work with current FCTD processing

                    % Get CTD data
                    FCTD.time=ctd.ctddnum; %Currently just ctdtime (seconds since powered on - will need to change this once we have gps data)
                    FCTD.pressure=ctd.P;
                    FCTD.temperature=ctd.T;
                    FCTD.conductivity=ctd.C;

                    % Get altimeter data
                    if ~isempty(alt) && isfield(alt,'alttime')
                        FCTD.altDist=interp1(alt.altdnum,alt.dst,ctd.ctddnum);
                    else
                        FCTD.altTime=nan(length(ctd.ctddnum),1);
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
                        ucontime=linspace(ctd.ctddnum(1),ctd.ctddnum(end),length(ctd.ctddnum)*20);
                        FCTD.uConductivity=reshape(interp1(epsi.epsidnum,double(epsi.s2_count),ucontime),20,[])';
                        clear ucontime
                    else
                        FCTD.uConductivity=nan(length(ctd.ctddnum),20);
                        disp(['No uConductivity data ' myASCIIfiles(i).name]);
                    end


                    % If we want the fluorometer to be outputted as well then
                    % that is saved in shear channel 1
                    % Currently saving in the same format as uCond but we may
                    % not need it to be so high resolution
                    if ~isempty(epsi) && isfield(epsi,'s1_volt')  && ~isempty(ctd)
                        fluortime=linspace(ctd.ctddnum(1),ctd.ctddnum(end),length(ctd.ctddnum)*20);
                        FCTD.fluorometer=reshape(interp1(epsi.epsidnum,epsi.s1_volt,fluortime),20,[])';
                        clear fluortime
                    else
                        FCTD.fluorometer=nan(length(ctd.ctddnum),20);
                        disp(['No fluorometer data ' myASCIIfiles(i).name]);
                    end

                    % Get VectorNav data

                    % Add vnav.vnavdnum data, interpolate to CTD
                    if ~isempty(vnav) && isfield(vnav,'vnavtime')
                        for ix=1:3
                            FCTD.compass(:,ix)=interp1(vnav.vnavdnum,vnav.compass(:,ix),ctd.ctddnum);
                            FCTD.gyro(:,ix)=interp1(vnav.vnavdnum,vnav.gyro(:,ix),ctd.ctddnum);
                            FCTD.acceleration(:,ix)=(interp1(vnav.vnavdnum,vnav.acceleration(:,ix),ctd.ctddnum))./9.81;
                        end
                    else
                        FCTD.gyro=nan(length(ctd.ctddnum),3);
                        FCTD.acceleration=nan(length(ctd.ctddnum),3);
                        FCTD.compass=nan(length(ctd.ctddnum),3);
                    end

                    % Add GPS data

                    if ~isempty(GPS) && isfield(GPS,'gpstime')
                        FCTD.GPS.longitude=interp1(GPS.gpstime,GPS.longitude,ctd.ctddnum);
                        FCTD.GPS.latitude=interp1(GPS.gpstime,GPS.latitude,ctd.ctddnum);
                    else
                        FCTD.GPS.longitude=nan(length(ctd.ctddnum),1);
                        FCTD.GPS.latitude=nan(length(ctd.ctddnum),1);
                    end


                    % Save FCTD mat files to the new FCTD mat directory FCTDmat
                    save([FCTDdir  base '.mat'],'FCTD');
                    FastCTD_UpdateMATFileTimeIndex(FCTDdir,base,FCTD);
                    fprintf(1,'%s: Wrote  %s%s\n\n',datestr(now,'YY.mm.dd HH:MM:SS'), FCTDdir,myFCTDMATfile.name);
                    %                 if doGrid
                    %                     FCTD_GridData = Epsi_GridData(FCTD);
                    %                     save([GridDir base '.mat'],'FCTD_GridData');
                    %                     Epsi_UpdateMATFileTimeIndex(GridDir,base,FCTD_GridData);
                    %                     fprintf(1,'%s: Wrote  %s%s.mat\n\n',datestr(now,'YY.mm.dd HH:MM:SS'),GridDir,base);
                    %                 end
                end %end Bethan's addition for FCTD data
            end %end if doFCTD



        catch err
            disp(['So... this is the error for retranlating file ' myASCIIfiles(i).name]);
            disp(err);
            for j = 1:length(err.stack)
                disp([num2str(j) ' ' err.stack(j).name ' ' num2str(err.stack(j).line)]);
            end
            error('There was an error. See stack above')
        end

        % If the files are new then a new MAT file will be created
    elseif isempty(myMATfile)
        fprintf(1,'Translating %s%s\n',RawDir,myASCIIfiles(i).name);
        try
            disp([RawDir myASCIIfiles(i).name]);

            [epsi,ctd,alt,act,vnav,GPS] = mod_som_read_epsi_files_vGPS([RawDir myASCIIfiles(i).name],Meta_Data);
            % Add fileNum

            % Display file size, time, pressure, altimeter
            try
                disp('~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ')
                disp(['FILE SIZE: ' num2str(myASCIIfiles(i).bytes)])
                disp(['TIME: ' datestr(ctd.ctddnum(end))])
                disp(['PRESSURE: ' num2str(ctd.P(end))])
            catch
            end
            try
                disp(['ALTIMETER: ' num2str(alt.dst(end))])
            catch
            end


            if ~isempty(epsi) && isfield(epsi,'epsitime')
                save([MatDir  base '.mat'],'epsi','ctd','alt','vnav');
                Epsi_UpdateMATFileTimeIndex(MatDir,base,epsi.epsitime);
                fprintf(1,'%s: Wrote  %s%s.mat\n\n',datestr(now,'YY.mm.dd HH:MM:SS'),MatDir,base);
            elseif isempty(epsi) &&  ~isempty(ctd) && isfield(ctd,'ctdtime') %For the case where the is no epsi data, but there is ctd data
                save([MatDir  base '.mat'],'epsi','ctd','alt','vnav');
                Epsi_UpdateMATFileTimeIndex(MatDir,base,ctd.ctdtime);
                fprintf(1,'%s: Wrote  %s%s.mat\n\n',datestr(now,'YY.mm.dd HH:MM:SS'),MatDir,base);
            end

            if doFCTD
                % Bethan's addition for FCTD data
                if ~isempty(ctd) && isfield(ctd,'ctdtime')
                    %Create new structure FCTD with data from alt, ctd and epsi renamed to work with current FCTD processing

                    % Get CTD data
                    FCTD.time=ctd.ctddnum;
                    FCTD.pressure=ctd.P;
                    FCTD.temperature=ctd.T;
                    FCTD.conductivity=ctd.C;

                    % Get altimeter data
                    if ~isempty(alt) && isfield(alt,'alttime')
                        FCTD.altDist=interp1(alt.altdnum,alt.dst,ctd.ctddnum);
                    else
                        FCTD.altTime=nan(length(ctd.ctddnum),1);
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
                        ucontime=linspace(ctd.ctddnum(1),ctd.ctddnum(end),length(ctd.ctddnum)*20);
                        FCTD.uConductivity=reshape(interp1(epsi.epsidnum,double(epsi.s2_count),ucontime),20,[])';
                        clear ucontime
                    else
                        FCTD.uConductivity=nan(length(ctd.ctddnum),20);
                        disp(['No uConductivity data ' myASCIIfiles(i).name]);
                    end


                    % If we want the fluorometer to be outputted as well then
                    % that is saved in shear channel 1
                    % Currently saving in the same format as uCond but we may
                    % not need it to be so high resolution
                    if ~isempty(epsi) && isfield(epsi,'s1_volt')  && ~isempty(ctd)
                        fluortime=linspace(ctd.ctddnum(1),ctd.ctddnum(end),length(ctd.ctddnum)*20);
                        FCTD.fluorometer=reshape(interp1(epsi.epsidnum,epsi.s1_volt,fluortime),20,[])';
                        clear fluortime
                    else
                        FCTD.fluorometer=nan(length(ctd.ctddnum),20);
                        disp(['No fluorometer data ' myASCIIfiles(i).name]);
                    end

                    % Get VectorNav data


                    %Bethan to do:

                    if ~isempty(vnav) && isfield(vnav,'vnavtime')
                        for ix=1:3
                            FCTD.compass(:,ix)=interp1(vnav.vnavdnum,vnav.compass(:,ix),ctd.ctddnum);
                            FCTD.gyro(:,ix)=interp1(vnav.vnavdnum,vnav.gyro(:,ix),ctd.ctddnum);
                            FCTD.acceleration(:,ix)=(interp1(vnav.vnavdnum,vnav.acceleration(:,ix),ctd.ctddnum))./9.81;
                        end
                    else
                        FCTD.gyro=nan(length(ctd.ctddnum),3);
                        FCTD.acceleration=nan(length(ctd.ctddnum),3);
                        FCTD.compass=nan(length(ctd.ctddnum),3);
                    end


                    % Add GPS data

                    if ~isempty(GPS) && isfield(GPS,'gpstime')
                        FCTD.GPS.longitude=interp1(GPS.gpstime,GPS.longitude,ctd.ctddnum);
                        FCTD.GPS.latitude=interp1(GPS.gpstime,GPS.latitude,ctd.ctddnum);
                    else
                        FCTD.GPS.longitude=nan(length(ctd.ctddnum),1);
                        FCTD.GPS.latitude=nan(length(ctd.ctddnum),1);
                    end


                    % Save FCTD mat files to the new FCTD mat directory FCTDmat
                    save([FCTDdir  base '.mat'],'FCTD');
                    FastCTD_UpdateMATFileTimeIndex(FCTDdir,base,FCTD);
                    fprintf(1,'%s: Wrote  %s%s\n\n',datestr(now,'YY.mm.dd HH:MM:SS'), FCTDdir,myFCTDMATfile.name);
                    %                 if doGrid
                    %                     FCTD_GridData = Epsi_GridData(FCTD);
                    %                     save([GridDir base '.mat'],'FCTD_GridData');
                    %                     Epsi_UpdateMATFileTimeIndex(GridDir,base,FCTD_GridData);
                    %                     fprintf(1,'%s: Wrote  %s%s.mat\n\n',datestr(now,'YY.mm.dd HH:MM:SS'),GridDir,base);
                    %                 end

                end %end Bethan's addition for FCTD data
            end %end if doFCTD

        catch err
            disp(['So... this is the error for tranlating file ' myASCIIfiles(i).name]);
            disp(err);
            for j = 1:length(err.stack)
                disp([num2str(j) ' ' err.stack(j).name ' ' num2str(err.stack(j).line)]);
            end
            error('There was an error. See stack above.')
        end
    end
    clear FCTD
end

clear epsi ctd alt act vnav;
end

%produce matlab indexing file for faster loading
function Epsi_UpdateMATFileTimeIndex(dirname,filename,datatime)
if exist([dirname '/Epsi_MATfile_TimeIndex.mat'],'file')
    load([dirname '/Epsi_MATfile_TimeIndex.mat']);
    ind = strncmp(filename,Epsi_MATfile_TimeIndex.filenames,length(filename));
    if sum(ind) ~= 1
        Epsi_MATfile_TimeIndex.filenames = [Epsi_MATfile_TimeIndex.filenames; {filename}];
        Epsi_MATfile_TimeIndex.timeStart = cat(1,Epsi_MATfile_TimeIndex.timeStart,datatime(1));
        Epsi_MATfile_TimeIndex.timeEnd = cat(1,Epsi_MATfile_TimeIndex.timeEnd,datatime(end));
    else
        Epsi_MATfile_TimeIndex.timeStart(ind) = datatime(1);
        Epsi_MATfile_TimeIndex.timeEnd(ind) = datatime(end);
    end
else
    Epsi_MATfile_TimeIndex.filenames = {filename};
    Epsi_MATfile_TimeIndex.timeStart = datatime(1);
    Epsi_MATfile_TimeIndex.timeEnd = datatime(end);
end
save([dirname '/Epsi_MATfile_TimeIndex.mat'],'Epsi_MATfile_TimeIndex');
end


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
