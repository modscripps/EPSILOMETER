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
% set the TimeIndex to a nonsense time such that all new time values will
% definitely be later.
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

% Sort indices of new data files by time
if ~isempty(deployment_MATfile_TimeIndex.timeEnd)
    TimeMin = nanmax(deployment_MATfile_TimeIndex.timeEnd(end));
else
    TimeMin = datenum(-1,1,1,0,0,0);
end
ind = find(Epsi_MATfile_TimeIndex.timeEnd > TimeMin);

[~,iSorted] = sort(Epsi_MATfile_TimeIndex.timeEnd(ind));
ind = ind(iSorted);

% Let's estimate how many new rows of data you'll need.

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
            if isstruct(epsi) && ~isempty(epsi.epsitime) && (epsi.epsitime(end)>=TimeMin) %now = Epsi_GUI_data.currentTime
                epsiOut = Epsi_MergeMatFiles(epsiOut,epsi);
                deployment_MATfile_TimeIndex.timeStart(end+1) = Epsi_MATfile_TimeIndex.timeStart(ind(i));
                deployment_MATfile_TimeIndex.timeEnd(end+1) =  Epsi_MATfile_TimeIndex.timeEnd(ind(i));
                deployment_MATfile_TimeIndex.filenames{end+1} = Epsi_MATfile_TimeIndex.filenames{ind(i)};
                
                clear epsi;
                
                epsi = epsiOut;
                
                % Make sure data are in order of increasing time and have unique time
                % values
                % 'unique' finds the unique values and returns them in sorted order.
                epsiFields = fields(epsi);
                [~,iSorted] = unique(epsi.epsitime);
                for iField=1:length(epsiFields)
                    epsi.(epsiFields{iField}) = epsi.(epsiFields{iField})(iSorted);
                end
                
                save(fullfile(Meta_Data.Epsipath,['deployment_MATfile_TimeIndex']),'deployment_MATfile_TimeIndex','-v7.3');
                save(fullfile(Meta_Data.Epsipath,['epsi_' Meta_Data.deployment '.mat']),'epsi','-v7.3');
                
            end
        end
        % Ctd
        % ---------------
        try
            load([Meta_Data.MATpath '/' Epsi_MATfile_TimeIndex.filenames{ind(i)} '.mat'],'ctd');
        catch err
            disp(err);
            disp([Meta_Data.MATpath '/' Epsi_MATfile_TimeIndex.filenames{ind(i)} '.mat']);
        end
        
        if exist('ctd','var') && ~isempty(ctd)
            if isstruct(ctd) && ~isempty(ctd.ctdtime) && (ctd.ctdtime(end)>=TimeMin) %now = Epsi_GUI_data.currentTime
                ctdOut = Epsi_MergeMatFiles(ctdOut,ctd);
            end
            clear ctd;
            
            
            ctd = ctdOut;
            
            % Make sure data are in order of increasing time and have unique time
            % values
            % 'unique' finds the unique values and returns them in sorted order.
            ctdFields = fields(ctd);
            [~,iSorted] = unique(ctd.ctdtime);
            for iField=1:length(ctdFields)
                ctd.(ctdFields{iField}) = ctd.(ctdFields{iField})(iSorted);
            end
            save(fullfile(Meta_Data.CTDpath, ['ctd_' Meta_Data.deployment '.mat']),'ctd','-v7.3');
        end
        % Alt 
        % ---------------
        try
            load([Meta_Data.MATpath '/' Epsi_MATfile_TimeIndex.filenames{ind(i)} '.mat'],'alt');
        catch err
            disp(err);
            disp([Meta_Data.MATpath '/' Epsi_MATfile_TimeIndex.filenames{ind(i)} '.mat']);
        end
        
        if exist('alt','var') && ~isempty(alt)
            if isstruct(alt) && ~isempty(alt.alttime) && (alt.alttime(end)>=TimeMin) %now = Epsi_GUI_data.currentTime
                altOut = Epsi_MergeMatFiles(altOut,alt);
            end
            clear alt;
        
        
        
        alt = altOut;
        
        % Make sure data are in order of increasing time and have unique time
        % values
        % 'unique' finds the unique values and returns them in sorted order.
        altFields = fields(alt);
        [~,iSorted] = unique(alt.alttime);
        for iField=1:length(altFields)
            alt.(altFields{iField}) = alt.(altFields{iField})(iSorted);
        end
        
        save(fullfile(Meta_Data.CTDpath, ['alt_' Meta_Data.deployment '.mat']),'alt','-v7.3');
        
        end
    end %end loop through ind
end %end if ~isempty(ind)










