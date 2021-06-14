function [epsiOut,ctdOut,altOut] = merge_new_mat_files(Meta_Data)

% Load time index of .ascii files that have been converted to .mat
try
    load([Meta_Data.MATpath '/' '/Epsi_MATfile_TimeIndex.mat']);
catch err
    if strcmp(err.identifier,'MATLAB:nonExistentField') || strcmp(err.identifier,'MATLAB:load:couldNotReadFile')
        error('There is no Epsi_MATfile_TimeIndex.mat')
    else
        error('Error loading Epsi_MATfile_TimeIndex.mat')
    end
end

% Load time index of .mat files that have been merged into
% epsi_depl*.mat,ctd_depl*.mat,alt_depl*.mat, etc. If it doesn't exist yet,
% set the TimeIndex to a nonsense time that all new data will be after
try
    load([Meta_Data.Epsipath '/deployment_MATfile_TimeIndex.mat']);
catch err
    if strcmp(err.identifier,'MATLAB:nonExistentField') || strcmp(err.identifier,'MATLAB:load:couldNotReadFile')
        deployment_MATfile_TimeIndex.timeStart = [];
        deployment_MATfile_TimeIndex.timeEnd = [];
        deployment_MATfile_TimeIndex.filenames = {};
    else
        error('Error loading deployment_MATfile_TimeIndex.mat')
    end
end

% Load epsi, ctd, and alt deployment data if it exists
if exist(fullfile(Meta_Data.Epsipath, ['epsi_' Meta_Data.deployment '.mat']),'file')
    data = load(fullfile(Meta_Data.Epsipath, ['epsi_' Meta_Data.deployment '.mat']));
    epsiOut = data.epsi;
else
    epsiOut = [];
end
if exist(fullfile(Meta_Data.CTDpath, ['ctd_' Meta_Data.deployment '.mat']),'file')
    data = load(fullfile(Meta_Data.CTDpath, ['ctd_' Meta_Data.deployment '.mat']));
    ctdOut = data.ctd;
else
    ctdOut = [];
end
if exist(fullfile(Meta_Data.CTDpath, ['alt_' Meta_Data.deployment '.mat']),'file')
    data = load(fullfile(Meta_Data.CTDpath, ['alt_' Meta_Data.deployment '.mat']));
    altOut = data.alt;
else
    altOut = [];
end

% Define indices of new data
if ~isempty(deployment_MATfile_TimeIndex.timeEnd)
    TimeMin = deployment_MATfile_TimeIndex.timeEnd(end);
else
    TimeMin = datenum(-1,1,1,0,0,0);
end
ind = find(Epsi_MATfile_TimeIndex.timeEnd > TimeMin);

[~,iSorted] = sort(Epsi_MATfile_TimeIndex.timeEnd(ind));
ind = ind(iSorted);


if ~isempty(ind)
    for i = 1:length(ind)
        
        % Epsi
        % ---------------
        try
            load([Meta_Data.MATpath '/' Epsi_MATfile_TimeIndex.filenames{ind(i)} '.mat'],'epsi');
        catch err
            disp(err);
            disp([Meta_Data.MATpath '/' Epsi_MATfile_TimeIndex.filenames{ind(i)} '.mat']);
        end
        
        if exist('epsi','var')
            if isstruct(epsi) && ~isempty(epsi.epsitime) && (epsi.epsitime(end)>=TimeMin && epsi.epsitime(1) <= now) %now = Epsi_GUI_data.currentTime
                epsiOut = Epsi_MergeMatFiles(epsiOut,epsi);
                deployment_MATfile_TimeIndex.timeStart(end+1) = Epsi_MATfile_TimeIndex.timeStart(ind(i));
                deployment_MATfile_TimeIndex.timeEnd(end+1) =  Epsi_MATfile_TimeIndex.timeEnd(ind(i));
                deployment_MATfile_TimeIndex.filenames{end+1} = Epsi_MATfile_TimeIndex.filenames{ind(i)};
            end
            clear epsi;
        end
        
        epsi = epsiOut;
        save(fullfile(Meta_Data.Epsipath,['deployment_MATfile_TimeIndex']),'deployment_MATfile_TimeIndex','-v7.3');
        save(fullfile(Meta_Data.Epsipath,['epsi_' Meta_Data.deployment '.mat']),'epsi','-v7.3');
        
        % Ctd
        % ---------------
        try
            load([Meta_Data.MATpath '/' Epsi_MATfile_TimeIndex.filenames{ind(i)} '.mat'],'ctd');
        catch err
            disp(err);
            disp([Meta_Data.MATpath '/' Epsi_MATfile_TimeIndex.filenames{ind(i)} '.mat']);
        end
        
        if exist('ctd','var')
            if isstruct(ctd) && ~isempty(ctd.ctdtime) && (ctd.ctdtime(end)>=TimeMin && ctd.ctdtime(1) <= now) %now = Epsi_GUI_data.currentTime
                ctdOut = Epsi_MergeMatFiles(ctdOut,ctd);
            end
            clear ctd;
        end
        
        ctd = ctdOut;
        save(fullfile(Meta_Data.CTDpath, ['ctd_' Meta_Data.deployment '.mat']),'ctd','-v7.3');
         
        % Alt
        % ---------------
        try
            load([Meta_Data.MATpath '/' Epsi_MATfile_TimeIndex.filenames{ind(i)} '.mat'],'alt');
        catch err
            disp(err);
            disp([Meta_Data.MATpath '/' Epsi_MATfile_TimeIndex.filenames{ind(i)} '.mat']);
        end
        
        if exist('alt','var')
            if isstruct(alt) && ~isempty(alt.alt.time) && (alti.altitime(end)>=TimeMin && alt.alttime(1) <= now) %now = Epsi_GUI_data.currentTime
                altOut = Epsi_MergeMatFiles(altOut,alt);
            end
            clear alt;
        end
        
        
        alt = altOut;
        save(fullfile(Meta_Data.CTDpath, ['alt_' Meta_Data.deployment '.mat']),'alt','-v7.3');
        
    end %end loop through ind
end %end if ~isempty(ind)



