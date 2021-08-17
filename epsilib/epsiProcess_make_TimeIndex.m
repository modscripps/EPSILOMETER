function [TimeIndex] = epsiProcess_make_TimeIndex(MATdir,str_to_match)
% [TimeIndex] = epsiProcess_make_TimeIndex(MATdir)
%
% Remake MATfile TimeIndex from scratch. Load all the files and
% find the first and last timestamps
%
% ex) MATdir = '/Volumes/MOD_data_1-1/FCTD_EPSI/RAW_0715/mat';
%     str_to_match = 'EPSI_B*.mat';
%     fileList = dir(fullfile(MATdir,str_to_match));

fileList = dir(fullfile(MATdir,str_to_match));
for iF=1:length(fileList)
    load(fullfile(MATdir,fileList(iF).name),'epsi');
    filename = fileList(iF).name(1:end-4);
    epsiProcess_update_TimeIndex(MATdir,filename,epsi);
end


