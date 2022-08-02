function obj = epsiProcess_gridProfiles(obj,z)

% NC 9/9/21 To do - This could be faster. Right now it always interpolates
% every profile and concatenates them. Could check to see if current
% profile already exists in grid

% Try loading griddedProfiles. If it doesn't exist, we'll
% make it
if exist(fullfile(obj.Meta_Data.paths.profiles,'griddedProfiles.mat'),'file')==2
    grid_exists = 1;
    G = load(fullfile(obj.Meta_Data.paths.profiles,'griddedProfiles'));
    grid = G.grid;
else
    grid_exists = 0;
end

fileList = dir(fullfile(obj.Meta_Data.paths.profiles,'Profile*.mat'));
for iFile=1:length(fileList)

    load(fullfile(obj.Meta_Data.paths.profiles,fileList(iFile).name))

    if isfield(Profile,'pr') && ~isempty(Profile.pr)%Only add to

        % Interpolate this profile to standard pressure grid
        gridNew = epsiProcess_interpolate_Profile_to_P(Profile,z);

        % Add gridded profile to full grid. If full grid doesn't
        % exist yet, create it.
        if grid_exists
            varList = fields(grid);
            for iVar=1:length(varList)
                if ~strcmp(varList{iVar},'mission') && ~strcmp(varList{iVar},'vehicle_name') && ~strcmp(varList{iVar},'deployment')
                    grid.(varList{iVar})(:,end+1) = gridNew.(varList{iVar});
                end

            end
        else
            varList = fields(gridNew);
            for iVar=1:length(varList)
                grid.(varList{iVar})(:,1) = gridNew.(varList{iVar});
            end
            grid_exists = 1;
        end

        % Keep only unique profiles - if you're running this in
        % realtime you could have part of a profiles followed by the
        % rest of it.
        % ... but keep the last one if there is more than one.
        % TO DO: make this step less janky
        profNumX = fliplr(1:length(grid.profNum));
        profNumY = fliplr(grid.profNum);
        [~,u] = unique(profNumY);
        for iVar=1:length(varList)
            if ~strcmp(varList{iVar},'mission') && ~strcmp(varList{iVar},'vehicle_name') ...
                    && ~strcmp(varList{iVar},'deployment') && ~strcmp(varList{iVar},'pr') ...
                    && ~strcmp(varList{iVar},'z')
                grid.(varList{iVar}) = grid.(varList{iVar})(:,profNumX(u));
            end
        end

        close

        grid.pr = grid.pr(:,1);
        grid.z = grid.z(:,1);
        grid.mission = grid.mission(:).';
        grid.vehicle_name = grid.vehicle_name(:).';
        grid.deployment = grid.deployment(:).';

        saveName = fullfile(obj.Meta_Data.paths.profiles,'griddedProfiles.mat');
        save(saveName,'grid');

    end
end