fileList = dir(fullfile(obj.Meta_Data.L1path,'Profile*.mat'));
for iFile=1:length(fileList)
    
    load(fullfile(obj.Meta_Data.L1path,fileList(iFile).name))
    % Interpolate this profile to standard pressure grid
    gridPnew = obj.f_interpolateProfileToP(Profile,P);
    
    % Try loading griddedProfiles. If it doesn't exist, we'll
    % make it
    if exist(fullfile(obj.Meta_Data.L1path,'griddedProfiles.mat'),'file')==2
        load(fullfile(obj.Meta_Data.L1path,'griddedProfiles'));
    end
    
    % Add gridded profile to full grid. If full grid doesn't
    % exist yet, create it.
    if exist('gridP','var')
        varList = fields(gridP);
        for iVar=1:length(varList)
            gridP.(varList{iVar})(:,end+1) = gridPnew.(varList{iVar});
            
        end
    else
        varList = fields(gridPnew);
        for iVar=1:length(varList)
            gridP.(varList{iVar})(:,1) = gridPnew.(varList{iVar});
        end
    end
    
    % Keep only unique profiles - if you're running this in
    % realtime ou could have part of a profiles followed by the
    % rest of it.
    % ... but keep the last one if there is more than one.
    % TO DO: make this step less janky
    profNumX = fliplr(1:length(gridP.profNum));
    profNumY = fliplr(gridP.profNum);
    [~,u] = unique(profNumY);
    for iVar=1:length(varList)
        gridP.(varList{iVar}) = gridP.(varList{iVar})(:,profNumX(u));
    end
    
    close
end

saveName = fullfile(obj.Meta_Data.L1path,'griddedProfiles');
eval(['save ' saveName ' gridP']);
