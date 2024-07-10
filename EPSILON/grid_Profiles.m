function [gridProfiles] = grid_Profiles(Meta_Data)
% [gridProfiles] = grid_Profiles(Meta_Data)
%
% Make 2-D arrays of profile data from all Profiles in deployment.
%
% Nicole Couto | August 2020
% -------------------------------------------------------------------------
%% Get list of profiles in the deployment
% -------------------------------------------------------------------------
profList = list_profiles_in_dir(Meta_Data.paths.profiles);

% Number of profiles in the deployment
nProfs = length(profList);

%% Load the first profile to get variable names
load(fullfile(Meta_Data.paths.profiles,profList{1}));
%profNum = str2double(profList{1}(8:10));
%eval(['Profile = ' sprintf('Profile%03.0f',profNum) ';']);

% Add average values of power spectral densities and coherences between 10-45 Hz
Profile = add_avg_psd_to_profile(Profile);


%% Make lists of variable names
% -------------------------------------------------------------------------
fieldNames = fields(Profile);

% Get rid of the fields that have the scan indices - no need to grid those
idxToss = contains(fieldNames,'ind_range');
fieldNames = fieldNames(~idxToss);

% Find the fields that are either 1- or 2-columns with length==nbscan
idxVars = false(length(fieldNames),1);
for iField=1:length(fieldNames)
    sz = size(Profile.(fieldNames{iField}));
    if any(sz==Profile.nbscan) &&  sum(sz>2)==1
        idxVars(iField) = true;
    end
end
variableList = fieldNames(idxVars);

% Some fields are inside structures. Loop through all the first level
% fields and find the second-level fields with 1- or 2-columns and
% length==nbscan
variableList2 = {};
structIdx = find(cell2mat(cellfun(@(F) isstruct(Profile.(F)),fieldNames,...
    'uniformoutput',0)));
for iStruct=1:length(structIdx)
    fieldNames2 = fields(Profile.(fieldNames{structIdx(iStruct)}));
    idxVars2 = false(length(fieldNames2),1);
    for iField=1:length(fieldNames2)
        sz = size(Profile.(fieldNames{structIdx(iStruct)}).(fieldNames2{iField}));
        if any(sz==Profile.nbscan) &&  sum(sz>2)==1
            idxVars2(iField) = true;
        end
    end

    if sum(idxVars2)>0
        structList = repmat({fieldNames{structIdx(iStruct)}},sum(idxVars2),1);
        variableList2 = [variableList2;...
            structList, fieldNames2(idxVars2)];
    end
end

%% Create pressure and time arrays
% -------------------------------------------------------------------------
% Since we don't yet know how long the longest profile is, make it
% unreasonably long at first. Later, we'll get rid of all the rows that are
% empty.
dz = Meta_Data.PROCESS.dz;
gridProfiles.pr(:,1) = 0:dz:1000;
nDepths = length(gridProfiles.pr);

% Allocate space to save the average dnum for each cast
gridProfiles.dnum = nan(1,nProfs);

%% Preallocate arrays to populate with data from each profile
% -------------------------------------------------------------------------
% Preallocate space for 1-column arrays, 2-column arrays, and cells the same
% length as nbscan
for iVar=1:length(variableList)
    currData = Profile.(variableList{iVar});

    % 1-column arrays
    if (isa(currData,'double') || isa(currData,'logical')) && any(size(currData)==1)

        switch variableList{iVar}
            case {'pr','dnum'}
                % do nothing for pr or dnum. You've already made those
                % arrays
            otherwise
                gridProfiles.(variableList{iVar}) = nan(nDepths,nProfs);
        end

        % 2-column arrays
    elseif (isa(currData,'double') || isa(currData,'logical')) && ~any(size(currData)==1)
        gridProfiles.([variableList{iVar},'_1']) = nan(nDepths,nProfs);
        gridProfiles.([variableList{iVar},'_2']) = nan(nDepths,nProfs);

        % cells
    elseif isa(currData,'cell')
        gridProfiles.(variableList{iVar}){nDepths,nProfs} = [];

    end
end %end loop through variableList

% Preallocate space for the 1-column arrays that have the correct size, but
% are inside structures
varLevel1 = unique(variableList2(:,1));
for iVar1=1:length(varLevel1)
    currLevel = contains(variableList2(:,1),varLevel1{iVar1});
    varLevel2 = variableList2(currLevel,2);
    for iVar2=1:length(varLevel2)
        gridProfiles.(varLevel1{iVar1}).(varLevel2{iVar2}) = nan(nDepths,nProfs);
    end
end


% Add a field for scan #
gridProfiles.iScan = nan(nDepths,nProfs);
variableList = [variableList; 'iScan'];

clear Profile*


%% Loop through all the profiles in the deployment and populate the arrays
% -------------------------------------------------------------------------
for iProf=1:length(profList)

    % Load profile
    load(fullfile(Meta_Data.paths.profiles,profList{iProf}));
    %profNum = str2double(profList{iProf}(8:10));
    %eval(['Profile = ' sprintf('Profile%03.0f',profNum) ';']);

    % Add average values of power spectral densities and coherences between 10-45 Hz
    try
    Profile = add_avg_psd_to_profile(Profile);
    catch
        disp('issue with profile')
    end


    % Add dummy iScan field so you can populate grid data later
    Profile.iScan = nan;

    % Find the indices of gridProfiles.pr that correspond to Profile.pr
    idx = ismember(gridProfiles.pr,Profile.pr);
    % Make sure pressure array in this profile is in increments of dz and
    % that pressure idices are sequential
    if ~(all(diff(Profile.pr)==dz) && all(diff(find(idx))==1))

        error('Profiles are not equally spaced. It might be time to rewrite gridProfiles to interpolate data to common pressure array.')

    elseif all(diff(Profile.pr)==dz) && all(diff(find(idx))==1)

        % Loop through variables not inside structures
        for iVar=1:length(variableList)

            currData = Profile.(variableList{iVar});

            % 1-column arrays
            if (isa(currData,'double') || isa(currData,'logical')) && any(size(currData)==1)
                switch variableList{iVar}
                    case 'pr'
                        % do nothing for pressure. That's a single common array
                    case 'dnum'
                        % Get the average dnum for the profile
                        gridProfiles.dnum(iProf) = nanmean(Profile.(variableList{iVar}));
                    case 'iScan'
                        gridProfiles.(variableList{iVar})(idx,iProf) = 1:Profile.nbscan;
                    otherwise
                        gridProfiles.(variableList{iVar})(idx,iProf) = Profile.(variableList{iVar});
                end

                % 2-column arrays
            elseif (isa(currData,'double') || isa(currData,'logical')) && ~any(size(currData)==1)
                gridProfiles.([variableList{iVar},'_1'])(idx,iProf) = Profile.(variableList{iVar})(:,1);
                gridProfiles.([variableList{iVar},'_2'])(idx,iProf) = Profile.(variableList{iVar})(:,2);

                % cells
            elseif isa(currData,'cell')
                gridProfiles.(variableList{iVar})(idx,iProf) = Profile.(variableList{iVar});

            end
        end %end loop through variableList


        % Loop through variables inside structures and keep the ones that
        % have one size element == 1
        for iVar=1:length(variableList2)
            varLevel1 = unique(variableList2(:,1));
            for iVar1=1:length(varLevel1)
                currLevel = contains(variableList2(:,1),varLevel1{iVar1});
                varLevel2 = variableList2(currLevel,2);
                for iVar2=1:length(varLevel2)
                    gridProfiles.(varLevel1{iVar1}).(varLevel2{iVar2})(idx,iProf) = Profile.(varLevel1{iVar1}).(varLevel2{iVar2});
                end
            end
        end %end loop through variableList2


    end %end check if pressure array for this profile is okay

    clear Profile*

end %end loop through profList

% Now get rid of all the excess rows
excessRows = all(isnan(gridProfiles.iScan),2);

gridVars = fields(gridProfiles);
for iVar=1:length(gridVars)
    switch gridVars{iVar}
        case 'dnum'
            % Do nothing for dnum
        case varLevel1 %for any of the strucures, loop through variable inside
            currLevel = contains(variableList2(:,1),gridVars{iVar});
            varLevel2 = variableList2(currLevel,2);
            for iVar2=1:length(varLevel2)
                gridProfiles.(gridVars{iVar}).(varLevel2{iVar2})(excessRows,:) = [];
            end
        otherwise
            gridProfiles.(gridVars{iVar})(excessRows,:) = [];
    end

end

% Add Meta_Data
gridProfiles.Meta_Data = Meta_Data;

%% Add varInfo
gridProfiles.varInfo.pr = {'CTD pressure','db'};
gridProfiles.varInfo.dnum = {'datenum','Matlab datenum'};
gridProfiles.varInfo.w = {'fall speed','db s^{-1}'};
gridProfiles.varInfo.t = {'temperature','C'};
gridProfiles.varInfo.s = {'salinity','psu'};
gridProfiles.varInfo.kvis = {'kinematic viscosity',''};
gridProfiles.varInfo.epsilon = {'turbulent kinetic energy dissipation rate calculated from Ps_shear_k', ''};
gridProfiles.varInfo.epsilon_co = {'turbulent kinetic energy dissipation rate calculated from Ps_shear_co_k', ''};
gridProfiles.varInfo.chi = {'temperature gradient dissipation rate','°C^2 s^{-1}'};
gridProfiles.varInfo.sh_fc = {'shear cutoff frequency, 1=uncorrected, 2=coherence-corrected', 'Hz'};
gridProfiles.varInfo.tg_fc = {'temperature gradient cutoff frequency, 1=uncorrected, 2=coherence-corrected', 'Hz'};
gridProfiles.varInfo.flag_tg_fc = {'temperature gradient cut off frequency is very high','0/1'};
gridProfiles.varInfo.Coh_10_45Hz_avg = {'average coherences between 10 and 45 Hz',''};
gridProfiles.varInfo.Pa_g_10_45Hz_avg = {'average acceleration power spectral densities between 10 and 45 Hz',''};
gridProfiles.varInfo.Ps_volt_10_45Hz_avg = {'average shear power spectral densities between 10 and 45 Hz',''};
gridProfiles.varInfo.Pt_volt_10_45Hz_avg = {'average temperature power spectral densities between 10 and 45 Hz',''};
gridProfiles.varInfo.iScan = {'scan indices',''};

%% Save data
save(fullfile(Meta_Data.paths.profiles,'gridded_Profiles'),'gridProfiles')
