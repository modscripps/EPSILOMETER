fileList = dir('Profile*.mat');
P  = 0:1:2200;
for iFile=1:length(fileList)
    
    load(fullfile(pwd,fileList(iFile).name));
                % Interpolate this profile to standard pressure grid
                gridPnew = interpolateProfileToP(Profile,P);


  
                    varList = fields(gridPnew);

                        for iVar=1:length(varList)
                            gridP.(varList{iVar})(:,iFile) = gridPnew.(varList{iVar});
                        end
                   

                
                % Keep only unique profiles - if you're running this in
                % realtime you could have part of a profile followed by the
                % rest of it.
                % ... but keep the last one if there is more than one.
                % TO DO: make this step less janky
                profNumX = fliplr(1:length(gridP.profNum));
                profNumY = fliplr(gridP.profNum);
                [~,u] = unique(profNumY);
                for iVar=1:length(varList)
                    gridP.(varList{iVar}) = gridP.(varList{iVar})(:,profNumX(u));
                end

end
            
save griddedProfiles gridP