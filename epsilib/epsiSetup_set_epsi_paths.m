function [Meta_Data] = set_epsi_paths(Meta_Data)

Meta_Data.datapath = pwd;

% NC - Stop making sd_raw directory. We don't use it, and more importantly
% it gets in the way of other functions searching the main directory for
% things called '*_raw'.
Meta_Data.RAWpath  = fullfile(Meta_Data.datapath,'raw');
Meta_Data.MATpath   = fullfile(Meta_Data.datapath,'mat');
% Meta_Data.CTDpath  = fullfile(Meta_Data.datapath,'merged');
% Meta_Data.Epsipath = fullfile(Meta_Data.datapath,'merged');
% Meta_Data.Mergedpath = fullfile(Meta_Data.datapath,'merged');
Meta_Data.L1path   = fullfile(Meta_Data.datapath,'profiles');
Meta_Data.Profilepath   = fullfile(Meta_Data.datapath,'profiles');
Meta_Data.FIGpath   = fullfile(Meta_Data.datapath,'figs');

% Create paths if they don't already exist
if ~exist(Meta_Data.L1path,'dir')
    eval([ '!mkdir ' strrep(Meta_Data.L1path,' ','\ ')]);
    disp(['Created new directory: ' Meta_Data.L1path])
end
% if ~exist(Meta_Data.Epsipath,'dir')
%     eval([ '!mkdir ' strrep(Meta_Data.Epsipath,' ','\ ')]);
%     disp(['Created new directory: ' Meta_Data.Epsipath])
% end
% if ~exist(Meta_Data.CTDpath,'dir')
%     eval([ '!mkdir ' strrep(Meta_Data.CTDpath,' ','\ ')]);
%     disp(['Created new directory: ' Meta_Data.CTDpath])
% end
if ~exist(Meta_Data.RAWpath,'dir')
    eval([ '!mkdir ' strrep(Meta_Data.RAWpath,' ','\ ')]);
    disp(['Created new directory: ' Meta_Data.RAWpath])
end
if ~exist(Meta_Data.MATpath,'dir')
    eval([ '!mkdir ' strrep(Meta_Data.MATpath,' ','\ ')]);
    disp(['Created new directory: ' Meta_Data.MATpath])
end
if ~exist(Meta_Data.FIGpath,'dir')
    eval([ '!mkdir ' strrep(Meta_Data.FIGpath,' ','\ ')]);
    disp(['Created new directory: ' Meta_Data.FIGpath])
end
% if ~exist(Meta_Data.Mergedpath,'dir')
%     eval([ '!mkdir ' strrep(Meta_Data.Mergedpath,' ','\ ')]);
%     disp(['Created new directory: ' Meta_Data.Mergedpath])
% end
if ~exist(Meta_Data.Profilepath,'dir')
    eval([ '!mkdir ' strrep(Meta_Data.Profilepath,' ','\ ')]);
    disp(['Created new directory: ' Meta_Data.Profilepath])
end