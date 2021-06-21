function Epsi_MakeMatFromRaw(dirs,Meta_Data,varargin)
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

% NC - make sure all dirs have a / at the end
for ii=1:numel(dirs)
    dirs{ii} = strrep([dirs{ii},'/'],'//','/');
%     % make windows compatible
%     dirs{ii} = strrep(dirs{ii}, '\', '/');
end

persistent argsNameToCheck;
if isempty(argsNameToCheck);
    argsNameToCheck = {'noSync','noGrid'};
end

rSync = true;
doGrid = false;

index = 1;
n_items = nargin-2;

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
    
myASCIIfiles = dir([RawDir, '*_raw']);
if (isempty(myASCIIfiles))
    myASCIIfiles = dir([RawDir, '*.ascii']);
end

for i=1:length(myASCIIfiles)
    base = myASCIIfiles(i).name(1:end-6);
    myMATfile = dir([MatDir base '.mat']);
    
    % if the MAT files are older than the data files, they will be retranslated
    if (~isempty(myMATfile) && datenum(myASCIIfiles(i).date)>datenum(myMATfile.date))
        fprintf(1,'Retranslating %s%s\n',MatDir,myMATfile.name);
        try
            disp([RawDir myASCIIfiles(i).name]);
            
            [epsi,ctd,alt] = mod_som_read_epsi_files_v2([RawDir myASCIIfiles(i).name],Meta_Data);
            if ~isempty(epsi) && isfield(epsi,'epsitime')
                save([MatDir  base '.mat'],'epsi','ctd','alt');
                Epsi_UpdateMATFileTimeIndex(MatDir,base,epsi);
                fprintf(1,'%s: Wrote  %s%s\n\n',datestr(now,'YY.mm.dd HH:MM:SS'), MatDir,myMATfile.name);
%                 if doGrid
%                     FCTD_GridData = Epsi_GridData(FCTD);
%                     save([GridDir base '.mat'],'FCTD_GridData');
%                     Epsi_UpdateMATFileTimeIndex(GridDir,base,FCTD_GridData);
%                     fprintf(1,'%s: Wrote  %s%s.mat\n\n',datestr(now,'YY.mm.dd HH:MM:SS'),GridDir,base);
%                 end
                
            end;
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
            
            [epsi,ctd,alt] = mod_som_read_epsi_files_v2([RawDir myASCIIfiles(i).name],Meta_Data);
            if ~isempty(epsi) && isfield(epsi,'epsitime')
                save([MatDir  base '.mat'],'epsi','ctd','alt');
                Epsi_UpdateMATFileTimeIndex(MatDir,base,epsi);
                fprintf(1,'%s: Wrote  %s%s.mat\n\n',datestr(now,'YY.mm.dd HH:MM:SS'),MatDir,base);
%                 if doGrid
%                     FCTD_GridData = Epsi_GridData(FCTD);
%                     save([GridDir base '.mat'],'FCTD_GridData');
%                     Epsi_UpdateMATFileTimeIndex(GridDir,base,FCTD_GridData);
%                     fprintf(1,'%s: Wrote  %s%s.mat\n\n',datestr(now,'YY.mm.dd HH:MM:SS'),GridDir,base);
%                 end
            end;
        catch err
            disp(['So... this is the error for tranlating file ' myASCIIfiles(i).name]);
            disp(err);
            for j = 1:length(err.stack)
                disp([num2str(j) ' ' err.stack(j).name ' ' num2str(err.stack(j).line)]);
            end
        end
    end;
end;

clear epsi ctd alt;
end

%produce matlab indexing file for faster loading
function Epsi_UpdateMATFileTimeIndex(dirname,filename,epsi)
if exist([dirname '/Epsi_MATfile_TimeIndex.mat'],'file')
    load([dirname '/Epsi_MATfile_TimeIndex.mat']);
    ind = strncmp(filename,Epsi_MATfile_TimeIndex.filenames,length(filename));
    if sum(ind) ~= 1
        Epsi_MATfile_TimeIndex.filenames = [Epsi_MATfile_TimeIndex.filenames; {filename}];
        Epsi_MATfile_TimeIndex.timeStart = cat(1,Epsi_MATfile_TimeIndex.timeStart,epsi.epsitime(1));
        Epsi_MATfile_TimeIndex.timeEnd = cat(1,Epsi_MATfile_TimeIndex.timeEnd,epsi.epsitime(end));
    else
        Epsi_MATfile_TimeIndex.timeStart(ind) = epsi.epsitime(1);
        Epsi_MATfile_TimeIndex.timeEnd(ind) = epsi.epsitime(end);
    end
else
    Epsi_MATfile_TimeIndex.filenames = {filename};
    Epsi_MATfile_TimeIndex.timeStart = epsi.epsitime(1);
    Epsi_MATfile_TimeIndex.timeEnd = epsi.epsitime(end);
end
save([dirname '/Epsi_MATfile_TimeIndex.mat'],'Epsi_MATfile_TimeIndex');
end
