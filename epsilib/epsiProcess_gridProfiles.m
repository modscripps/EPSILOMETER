function obj = epsiProcess_gridProfiles(obj,z)

fileList = dir(fullfile(obj.Meta_Data.paths.profiles,'Profile*.mat'));
for iFile=1:length(fileList)
    
    load(fullfile(obj.Meta_Data.paths.profiles,fileList(iFile).name))
    % Interpolate this profile to standard pressure grid
    gridNew = obj.f_interpolateProfileToZ(Profile,z);
    
    % Try loading griddedProfiles. If it doesn't exist, we'll
    % make it
    if exist(fullfile(obj.Meta_Data.paths.profiles,'griddedProfiles.mat'),'file')==2
        load(fullfile(obj.Meta_Data.paths.profiles,'griddedProfiles'));
    end
    
    % Add gridded profile to full grid. If full grid doesn't
    % exist yet, create it.
    if exist('grid','var')
        varList = fields(grid);
        for iVar=1:length(varList)
            grid.(varList{iVar})(:,end+1) = gridNew.(varList{iVar});
            
        end
    else
        varList = fields(gridNew);
        for iVar=1:length(varList)
            grid.(varList{iVar})(:,1) = gridNew.(varList{iVar});
        end
    end
    
    % Keep only unique profiles - if you're running this in
    % realtime ou could have part of a profiles followed by the
    % rest of it.
    % ... but keep the last one if there is more than one.
    % TO DO: make this step less janky
    profNumX = fliplr(1:length(grid.profNum));
    profNumY = fliplr(grid.profNum);
    [~,u] = unique(profNumY);
    for iVar=1:length(varList)
        grid.(varList{iVar}) = grid.(varList{iVar})(:,profNumX(u));
    end
    
    close
end

saveName = fullfile(obj.Meta_Data.paths.profiles,'griddedProfiles');
eval(['save ' saveName ' grid']);

end