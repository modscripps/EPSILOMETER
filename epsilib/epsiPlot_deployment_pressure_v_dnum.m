function [CTD,ax] = epsiPlot_deployment_pressure_v_dnum(matDir)
%
% Merge all ctd pressure and dnum data from files in this directory and
% plot. This will allow you to pick out filenames of interest. Filenames
% inclue date and time in the title. 
%
% INPUTS
%   matDir - directory containing .mat files of with ctd structures
% -------------------------------------------------------------------------

fileList = dir(fullfile(matDir, '*.mat'));
nFiles = length(fileList);

CTD.dnum = nan(5000*nFiles,1);
CTD.ctdtime = nan(5000*nFiles,1);
CTD.P = nan(5000*nFiles,1);
for iFile=1:nFiles
    load(fullfile(matDir,fileList(iFile).name),'ctd');   

    CTD.dnum = [CTD.dnum;ctd.dnum];
    CTD.ctdtime = [CTD.ctdtime;ctd.ctdtime];
    CTD.P = [CTD.P;ctd.P];    
end

figure
gap = [0.1 0.05];
margV = [0.1 0.05];
margH = [0.08 0.05];
ax(1) = subtightplot(2,1,1,gap,margV,margH);
plot(CTD.dnum,CTD.P)
ax(1).YDir = 'reverse';
ax(1).YLabel.String = 'Pressure';
datetick(ax(1),'x','HH:MM')

ax(2) = subtightplot(2,1,2,gap,margV,margH);
plot(CTD.ctdtime,CTD.P)
ax(2).YDir = 'reverse';
ax(2).YLabel.String = 'Pressure';
ax(2).XLabel.String = 'ctdtime';

figureStamp(getFilename)

    
    
    
    

