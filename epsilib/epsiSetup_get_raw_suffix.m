function [Meta_Data] = epsiSetup_get_raw_suffix(Meta_Data)

% Determine the suffix of filenames in the raw directory.
% epsiProcess_convert_new_raw_to_mat requires this information to know which files to
% convert.

% List possible suffix options. They need to be actual suffixes. Nothing
% can exist after these characters because Epsi_MakeMatFromRaw creates
% files based on the name BEFORE the suffix.
suffixOptions = {'.ascii','_raw','.raw'};

% List files in raw directory
% rawDirContents = dir(Meta_Data.RAWpath);
%ALB change to check on datapath because I do not see how the processing
%can cp files in the RAW folder prior to that step
rawDirContents = dir(Meta_Data.datapath);

% Todo: If there is nothing in the raw directory, list files in the main data
% directory

% Count files that match each of the suffix options
for iOpt=1:length(suffixOptions)
    matchCell = strfind({rawDirContents.name}, suffixOptions{iOpt});
    matchArray = cellfun(@(C) ~isempty(C),matchCell);
    counts(iOpt) = sum(matchArray);
end

% We use the suffix option with the most match counts
[~,idxChoice] = max(counts);
Meta_Data.rawfileSuffix = suffixOptions{idxChoice};


%ALB I think we should copy raw file now in the raw folder
copyfile(fullfile(Meta_Data.datapath,['*' Meta_Data.rawfileSuffix]),Meta_Data.RAWpath)