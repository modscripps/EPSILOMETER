function FastCTD_MakeMatFromRaw_MISOBOB2019(dirs,varargin)
% FastCTD_MakeMatFromRaw(DIRS)
% Script to make FCTD "raw" (ASCII) matfiles
% Run: FastCTD_MakeMatFromRaw(DIRS)
% where DIRS should be
% DIRS = {RawDir;
%         RawDirAway;
%         MatDir;
%         GridDir;}
%
% Written by Jody Klymak
% Updated 2011 06 21 by San Nguyen
% Updated 2012 09 29 by San Nguyen for EquatorMix2012



persistent argsNameToCheck;
if isempty(argsNameToCheck);
    argsNameToCheck = {'noSync','noGrid'};
end

rSync = true;
doGrid = true;

index = 1;
n_items = nargin-1;

while (n_items > 0)
    argsMatch = strcmpi(varargin{index},argsNameToCheck);
    i = find(argsMatch,1);
    if isempty(i)
        error('MATLAB:FastCTD_MakeMatFromRaw:wrongOption','Incorrect option specified: %s', varargin{index});
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
    end
end

if rSync
    RawDir = [dirs{1}];
    RawDirDuplicate = [dirs{2}];
    MatDir = [dirs{3}];
    if doGrid
        GridDir = [dirs{4}];
    end
    
    % check for valid directories
    if ~exist(RawDirDuplicate,'dir');
        error('Cannot find remote dir: %s',RawDirDuplicate);
    end
else
    RawDir = [dirs{1}];
    MatDir = [dirs{2}];
    if doGrid
        GridDir = [dirs{3}];
    end
end

if ~exist(RawDir,'dir');
    error('Cannot find local RawDir: %s',RawDir);
end

if ~exist(MatDir,'dir');
    error('Cannot find local MatDir: %s',MatDir);
end

if doGrid && ~exist(GridDir,'dir');
    error('Cannot find local GridDir: %s',GridDir);
end


% rsync remote and local directories.  You must connect to the server
% to make this work.
if rSync
    com = sprintf('/usr/bin/rsync -av %s %s',RawDir,RawDirDuplicate);
    fprintf(1,'Running: %s\n',com);
    unix(com);
    fprintf(1,'Done\n');
    RawDir = RawDirDuplicate;
end

% we apply different calibration for SBE SN 0058 because the 2019
% calibration from sea bird causes temperature to have an offset
% SERIALNO=0058
% TCALDATE=10-dec-17
ctd_cal.ta0 = 8.646787e-04;
ctd_cal.ta1 = 2.862445e-04;
ctd_cal.ta2 = -2.466554e-06;
ctd_cal.ta3 = 2.121614e-07;
ctd_cal.cg = -9.727914e-01;
ctd_cal.ch = 1.412282e-01;
ctd_cal.ci = -2.100128e-04;
ctd_cal.cj = 3.888261e-05;
ctd_cal.pa0 = -9.985992e-02;
ctd_cal.pa1 = 4.516155e-03;
ctd_cal.pa2 = -1.856201e-11;
ctd_cal.ptempa0 = -5.789793e+01;
ctd_cal.ptempa1 = 5.622564e+01;
ctd_cal.ptempa2 = -9.682869e-01;
ctd_cal.ptca0 = 5.212094e+05;
ctd_cal.ptca1 = -2.710947e+00;
ctd_cal.ptca2 = 6.387558e-02;
ctd_cal.ptcb0 = 2.444488e+01;
ctd_cal.ptcb1 = -6.250000e-04;
ctd_cal.ptcb2 = 0.000000e+00;


ctd_cal.cpcor = -9.570000e-08;
ctd_cal.ctcor = 3.250000e-06;
ctd_cal.cslope = 1.000000e+00;

myASCIIfiles = dir([RawDir '*.ascii']);

for i=1:length(myASCIIfiles)
    base = myASCIIfiles(i).name(1:end-6);
    myMATfile = dir([MatDir base '.mat']);
    
    % if the MAT files are older than the data files, they will be retranslated
    if (~isempty(myMATfile) && datenum(myASCIIfiles(i).date)>datenum(myMATfile.date))
        fprintf(1,'Retranslating %s%s\n',MatDir,myMATfile.name);
        try
            disp([RawDir myASCIIfiles(i).name]);
            FCTD = FastCTD_ReadASCII_with_cal([RawDir myASCIIfiles(i).name],ctd_cal);
            if ~isempty(FCTD) && isfield(FCTD,'pressure')
                save([MatDir  base '.mat'],'FCTD');
                FastCTD_UpdateMATFileTimeIndex(MatDir,base,FCTD);
                fprintf(1,'%s: Wrote  %s%s\n\n',datestr(now,'YY.mm.dd HH:MM:SS'), MatDir,myMATfile.name);
                if doGrid
                    FCTD_GridData = FastCTD_GridData(FCTD);
                    save([GridDir base '.mat'],'FCTD_GridData');
                    FastCTD_UpdateMATFileTimeIndex(GridDir,base,FCTD_GridData);
                    fprintf(1,'%s: Wrote  %s%s.mat\n\n',datestr(now,'YY.mm.dd HH:MM:SS'),GridDir,base);
                end
                
            end
        catch err
            disp(['So... this is the error for retranlating file ' myASCIIfiles(i).name]);
            disp(err);
            for j = 1:length(err.stack)
                disp([num2str(j) ' ' err.stack(j).name ' ' num2str(err.stack(j).line)]);
            end
        end
    % the files are new then a few MAT file will be created
    elseif isempty(myMATfile)
        fprintf(1,'Translating %s%s\n',RawDir,myASCIIfiles(i).name);
        try
            disp([RawDir myASCIIfiles(i).name]);
            FCTD = FastCTD_ReadASCII_with_cal([RawDir myASCIIfiles(i).name],ctd_cal);
            if ~isempty(FCTD) && isfield(FCTD,'pressure')
                save([MatDir base '.mat'],'FCTD');
                FastCTD_UpdateMATFileTimeIndex(MatDir,base,FCTD);
                fprintf(1,'%s: Wrote  %s%s.mat\n\n',datestr(now,'YY.mm.dd HH:MM:SS'),MatDir,base);
                if doGrid
                    FCTD_GridData = FastCTD_GridData(FCTD);
                    save([GridDir base '.mat'],'FCTD_GridData');
                    FastCTD_UpdateMATFileTimeIndex(GridDir,base,FCTD_GridData);
                    fprintf(1,'%s: Wrote  %s%s.mat\n\n',datestr(now,'YY.mm.dd HH:MM:SS'),GridDir,base);
                end
            end
        catch err
            disp(['So... this is the error for tranlating file ' myASCIIfiles(i).name]);
            disp(err);
            for j = 1:length(err.stack)
                disp([num2str(j) ' ' err.stack(j).name ' ' num2str(err.stack(j).line)]);
            end
        end
    end
end

clear FCTD FCTD_GridData;
end

%produce matlab indexing file for faster loading
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
