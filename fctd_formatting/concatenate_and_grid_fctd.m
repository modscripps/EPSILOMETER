function [FCTDall,FCTDgrid] = concatenate_and_grid_fctd(fctd_mat_dir)

% Get FCTD files in current directory
f = dir(fullfile(fctd_mat_dir,'EPSI*.mat'));
file_list = {f.name};

% Make a list of datenums for each file
for iF=1:length(file_list)
    yyyy = str2num(['20' file_list{iF}(5:6)]);
    mm = str2num(file_list{iF}(8:9));
    dd = str2num(file_list{iF}(11:12));
    HH = str2num(file_list{iF}(14:15));
    MM = str2num(file_list{iF}(16:17));
    SS = str2num(file_list{iF}(18:19));
    file_dnums(iF) = datenum(yyyy,mm,dd,HH,MM,SS);
end

% Look for FCTDall so you can just add the new ones
if exist(fullfile(fctd_mat_dir,'FCTDall.mat'),'file')==2 %If FCTDall already exists
    % Load FCTDall
    load(fullfile(fctd_mat_dir,'FCTDall'))
    
    % What's the last datenum?
    last_dnum = FCTDall.time(end);
    
    % Keep any file newer than 10 minutes before last_dnum
    idx = file_dnums > last_dnum - days(minutes(10));
    file_list = file_list(idx);
    
else % If FCTDall doesn't already exist
    % Load the first file and add microconductivity data
    load(fullfile(fctd_mat_dir,file_list{1}),'FCTD');

    % %ALB hack TFO2024 to remove chi field and reprocess
    % FCTD=rmfield(FCTD,'chi');
    % FCTD=rmfield(FCTD,'chi2');

    FCTD = add_microconductivity(FCTD);
    
    % Save data with chi
    save(fullfile(fctd_mat_dir,file_list{1}),'FCTD');
    
    FCTDall = FCTD;
end

% Load the rest of the files and merge into FCTD all
for iF=2:length(file_list)
    fprintf(sprintf('Adding file %3.0f of %3.0f to FCTDall \n',iF,length(file_list)));
    
    % Load file and add microconductivity data
    load(fullfile(fctd_mat_dir,file_list{iF}),'FCTD');
    % %ALB hack TFO2024 to remove chi field and reprocess
    % FCTD=rmfield(FCTD,'chi');
    % FCTD=rmfield(FCTD,'chi2');

    FCTD = add_microconductivity(FCTD);
    
    % Save data with chi
    save(fullfile(fctd_mat_dir,file_list{iF}),'FCTD');
    
    % Concatenate new file into FCTDall
    FCTDall = FastCTD_MergeFCTD(FCTDall,FCTD);
end

% You might have reprocessed some of the same data. Find unique dnums and
% only keep those data
[~,iU] = unique(FCTDall.time);

% Get field list
field_names = fields(FCTDall);
for f=1:length(field_names)
    FCTDall.(field_names{f}) = FCTDall.(field_names{f})(iU,:);
end

% Save concatenated file
save(fullfile(fctd_mat_dir,'FCTDall'),'FCTDall')

%% Grid data and save
% NC May 2023 - grid both up and downcast data
% Get downcasts and upcasts

FCTDdown = FastCTD_GridData(FCTDall,'downcast');
FCTDup = FastCTD_GridData(FCTDall,'upcast');

% Merge all the structure fields for the down and up casts
field_list = fields(FCTDdown);
for iF=1:length(field_list)
    % skip tGrid but do all the others
    if ~strcmp(field_list{iF},'tGrid') 
        FCTDboth.(field_list{iF}) = [FCTDdown.(field_list{iF}),...
                                     FCTDup.(field_list{iF})];
    end
end

% Now sort the merged profiles by dnum
[~,iSort] = sort(FCTDboth.time);
field_list = fields(FCTDboth);
for iF=1:length(field_list)
    % Don't sort depth by dnum, but do make it a Nx1 field
    if strcmp(field_list{iF},'depth')
        FCTDgrid.(field_list{iF}) = FCTDboth.(field_list{iF})(:,1);
    else
    FCTDgrid.(field_list{iF}) = FCTDboth.(field_list{iF})(:,iSort);
    end
end


% Save grid
save(fullfile(fctd_mat_dir,'FCTDgrid'),'FCTDgrid','FCTDup','FCTDdown')

end %end concatenate_fctd_data

% -------------------------------------------------------------------------function [FCTD] = add_microconductivity(FCTD)

function [FCTD] = add_microconductivity(FCTD)
% Check for chi data
if isfield(FCTD,'uConductivity') && ~isfield(FCTD,'chi')
    
    % Load chi_param
    chi_param=FCTD_DefaultChiParam;
    chi_param.fs=320;
    chi_param.min_spd=0.01; %TFO RR2410 upcast have a slow last part of up cast.
    
    % Reshape the data into a long vector
    FCTD.ucon=reshape(FCTD.uConductivity',20*length(FCTD.time),1)/2^24;
    
    % Make the higher resolution time vector for microconductivity
    FCTD.microtime=linspace(FCTD.time(1),FCTD.time(end),chi_param.fs./chi_param.fslow*length(FCTD.time))';
    
    % Convert microconductivity to chi
    FCTD = add_chi_microMHA_v2(FCTD,chi_param);

    
end

% Remove descriptive fields
if isfield(FCTD,'chi_param')
    FCTD = rmfield(FCTD,'chi_param');
end
if isfield(FCTD,'chi_meta')
    FCTD = rmfield(FCTD,'chi_meta');
end

end % end add_microconductivity
