function [PressureTimeseries] = epsiProcess_make_PressureTimeseries(MATdir,str_to_match)
% [PressureTimeseries] = epsiProcess_make_PressureTimeseries(MATdir,str_to_match)
%
% Remake pressure timeseries from scratch. Load all the files and
% concatenate the pressure data
%
% ex) MATdir = '/Volumes/MOD_data_1-1/FCTD_EPSI/RAW_0715/mat';
%     str_to_match = 'EPSI_B*.mat';
%     fileList = dir(fullfile(MATdir,str_to_match));

fileList = dir(fullfile(MATdir,str_to_match));
for iF=1:length(fileList)
    load(fullfile(MATdir,fileList(iF).name),'ctd');
    epsiProcess_update_PressureTimeseries(MATdir,ctd);
end


