function [profList] = list_profiles_in_dir(dataDir)
% [profList] = list_profiles_in_dir(dataDir)
%
% OUTPUTS
%   profList - a cell string of Profile###.mat in dataDir

dirContents = dir(dataDir);
fileNames = {dirContents(:).name};
idxProf = cell2mat(cellfun(@(C) ~isempty(regexp(C,'Profile[0-9][0-9][0-9].mat', 'once')),...
    fileNames,'UniformOutput',false));
profList = fileNames(idxProf);