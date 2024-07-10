function obj = epsiProcess_gridProfiles(obj,z)

% NC 9/9/21 To do - This could be faster. Right now it always interpolates
% every profile and concatenates them. Could check to see if current
% profile already exists in grid

% Try loading griddedProfiles. If it doesn't exist, we'll
% make it
if exist(fullfile(obj.Meta_Data.paths.profiles,'griddedProfiles.mat'),'file')==2
    grid_exists = 1;
    G = load(fullfile(obj.Meta_Data.paths.profiles,'griddedProfiles'));
    GRID = G.GRID;
else
    grid_exists = 0;
end
% Actually, always re-grid from the beginning since you might be
% changing the depth array
grid_exists = 0;
clear GRID

fileList = dir(fullfile(obj.Meta_Data.paths.profiles,'Profile*.mat'));
for iFile=1:length(fileList)

    disp(sprintf('Gridding %03.0f of %03.f',iFile,length(fileList)))
    load(fullfile(obj.Meta_Data.paths.profiles,fileList(iFile).name));

    if isfield(Profile,'pr') && ~isempty(Profile.pr)%Only add to

        % Interpolate this profile to standard pressure grid
        gridNew = epsiProcess_interpolate_Profile_to_P(Profile,z);

        % Add gridded profile to full GRID. If full grid doesn't
        % exist yet, create it.
        if grid_exists
            varList = fields(GRID);
            for iVar=1:length(varList)
                % temporarily ignores lat and lon to fix gps error
                if sum(strcmp(varList{iVar},{'mission','vehicle_name','deployment','filenames'}))==0
                    if (isfield(gridNew,varList{iVar}))
                        GRID.(varList{iVar})(:,end+1) = gridNew.(varList{iVar});
                    else
                        if(sum(strcmp(varList{iVar},{'bottom_depth','profNum','dnum','latitude','longitude'})))
                            GRID.(varList{iVar})(:,end+1) = NaN;
                        else
                            GRID.(varList{iVar})(:,end+1) = NaN(size(gridNew.z));
                        end
                    end
                end

            end
        else
            varList = fields(gridNew);
            for iVar=1:length(varList)
                GRID.(varList{iVar})(:,1) = gridNew.(varList{iVar});
            end
            grid_exists = 1;
        end

        % Keep only unique profiles - if you're running this in
        % realtime you could have part of a profiles followed by the
        % rest of it.
        % ... but keep the last one if there is more than one.
        % TO DO: make this step less janky
        profNumX = fliplr(1:length(GRID.profNum));
        profNumY = fliplr(GRID.profNum);
        [~,u] = unique(profNumY);
        for iVar=1:length(varList)
            if ~strcmp(varList{iVar},'mission') && ~strcmp(varList{iVar},'vehicle_name') ...
                    && ~strcmp(varList{iVar},'deployment') && ~strcmp(varList{iVar},'pr') ...
                    && ~strcmp(varList{iVar},'z') && ~strcmp(varList{iVar},'filenames')
                GRID.(varList{iVar}) = GRID.(varList{iVar})(:,profNumX(u));
            end
        end

        close

        GRID.pr = GRID.pr(:,1);
        GRID.z = GRID.z(:,1);
        GRID.mission = GRID.mission(:).';
        GRID.vehicle_name = GRID.vehicle_name(:).';
        GRID.deployment = GRID.deployment(:).';
        GRID.filenames{iFile} = Profile.filenames;
        
    saveName = fullfile(obj.Meta_Data.paths.profiles,'griddedProfiles.mat');
    save(saveName,'GRID');

    end %End if there is Profile.pr field

end %End loop through profiles