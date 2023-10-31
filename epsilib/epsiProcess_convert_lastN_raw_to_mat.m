function [matData] = epsiProcess_convert_lastN_raw_to_mat(RawDir,Meta_Data,N)
% epsiProcess_convert_lastN_raw_to_mat
%   - converts the last 
%
% Nicole Couto | Summer 2021
% --------------------------------------------------------
% INPUTS:
%   dirs = {RawDir;
%           RawDirAway;
%           MatDir;
%           FctdMatDir;}
%   Meta_Data: epsi Meta_Data structure
%   N = number of files to load and concatenate

if nargin<3
    N = 1;
end

% NC - Make matData for output even if there is no new data
matData.epsi = [];
matData.ctd = [];
matData.alt = [];
matData.act = [];
matData.vnav = [];
matData.gps = [];
matData.seg = [];
matData.spec = [];
matData.avgspec = [];
matData.dissrate = [];
matData.apf = [];
matData.fluor = [];
matData.ttv = [];


% Find files in RawDir ending in suffixSearch
suffixStr = Meta_Data.PROCESS.rawfileSuffix; %ex. *.raw, *.ascii, etc
suffixSearch = ['*' suffixStr];
myASCIIfiles = dir(fullfile(RawDir, suffixSearch));

% Loop through the last N files
for i=length(myASCIIfiles)-(N-1):length(myASCIIfiles)
    
    % Convert raw data to mat
    newData = mod_som_read_epsi_files_v4(fullfile(RawDir,myASCIIfiles(i).name),Meta_Data);
    use newData
    
    % Find data fields
    field_names = fields(newData);
    
    % Concatenate the data
    for iF=1:length(field_names)
        matData.(field_names{iF}) = epsiProcess_merge_mat_files(matData.(field_names{iF}),newData.(field_names{iF}));
    end
end

% Clean up matData - if a field is empty, get rid of it.
empty = structfun(@isempty, matData);
matData = rmfield(matData,field_names(empty));

end
