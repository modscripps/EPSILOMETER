function [Meta_Data] = set_epsi_paths(Meta_Data)

% ALB we moved the call earlier
% Meta_Data.paths.data = pwd;

% NC - Stop making sd_raw directory. We don't use it, and more importantly
% it gets in the way of other functions searching the main directory for
% things called '*_raw'.
Meta_Data.paths.raw_data  = fullfile(Meta_Data.paths.data,'raw');
Meta_Data.paths.mat_data   = fullfile(Meta_Data.paths.data,'mat');
Meta_Data.paths.profiles   = fullfile(Meta_Data.paths.data,'profiles');
Meta_Data.paths.figures   = fullfile(Meta_Data.paths.data,'figs');

suffixOptions = {'.ascii','_raw','.raw','.epsi','.bin','.modraw'};

% Count files that match each of the suffix options
for iOpt=1:length(suffixOptions)
    try
        copyfile(['*' suffixOptions{iOpt}],Meta_Data.paths.raw_data)
    catch
    end
end

% Create paths if they don't already exist
pathList = {'raw_data','mat_data','profiles','figures'};
for p=1:length(pathList)
    if ~exist(Meta_Data.paths.(pathList{p}),'dir')
        mkdir(Meta_Data.paths.(pathList{p}));
        %eval([ '!mkdir ' strrep(Meta_Data.paths.(pathList{p}),' ','\ ')]);
        disp(['Created new directory: ' Meta_Data.paths.(pathList{p})])
    end
end
