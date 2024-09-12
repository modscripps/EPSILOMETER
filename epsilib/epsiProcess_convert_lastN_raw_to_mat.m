function [matData] = epsiProcess_convert_lastN_raw_to_mat(RawDir,Meta_Data,N)
% epsiProcess_convert_lastN_raw_to_mat
%   - converts the last 
%
% aleboyer@ucsd.edu June2023 adding ttv and fluor
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
matData.epsi  = [];
matData.ctd   = [];
matData.alt   = [];
matData.act   = [];
matData.vnav  = [];
matData.gps   = [];
matData.fluor = [];
matData.ttv   = [];

% Find files in RawDir ending in suffixSearch
suffixStr = Meta_Data.PROCESS.rawfileSuffix; %ex. *.raw, *.ascii, etc
suffixSearch = ['*' suffixStr];
myASCIIfiles = dir(fullfile(RawDir, suffixSearch));

% Loop through the last N files
for i=length(myASCIIfiles)-(N-1):length(myASCIIfiles)
    
    % Convert raw data to mat
    newData = mod_som_read_epsi_files_v4(fullfile(RawDir,myASCIIfiles(i).name),Meta_Data);
    use newData
    
    % If you loaded more than one file, concatenate them. Otherwise, just
    % put into matData structure
    matData.epsi  = epsiProcess_merge_mat_files(matData.epsi,epsi);
    matData.ctd   = epsiProcess_merge_mat_files(matData.ctd,ctd);
    matData.alt   = epsiProcess_merge_mat_files(matData.alt,alt);
    matData.act   = epsiProcess_merge_mat_files(matData.act,act);
    matData.vnav  = epsiProcess_merge_mat_files(matData.vnav,vnav);
    matData.gps   = epsiProcess_merge_mat_files(matData.gps,gps);
    matData.ttv   = epsiProcess_merge_mat_files(matData.ttv,ttv);
    matData.fluor = epsiProcess_merge_mat_files(matData.fluor,fluor);
end

clear epsi ctd alt act vnav ttv fluor;
end
