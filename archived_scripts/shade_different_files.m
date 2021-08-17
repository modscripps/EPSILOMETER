function [] = shade_different_files(obj,timeOrSample)

% Add shaded boxes to a plot to denote different files
%
% USAGE:
%   [] = shade_different_files(epsi,timeOrSamplenum)
%
% INPUTS:
%   obj             = epsi_class object
%   timeOrSamplenum = 'time' if x-axis is time, or 'sample' if x-axis is
%                      sample number

%% Check for enough inputs
if nargin<2
    error('Choose "time" or "samplenum" for x-axis --> shade_different_files(dataStruct,timeOrSamplenum)')
end

%% Find where filenum changes
switch timeOrSample
    case {'epsitime','epsisample'}
        dataStruct = obj.epsi;
    case {'ctdtime','ctdsample'}
        dataStruct = obj.ctd;
end

dF = find(diff(dataStruct.fileNum)~=0);
firsts = [1;dF];
lasts = [dF+1;length(dataStruct.fileNum)];

%% Get filenames
fileDir = obj.filename;
fileList = dir(fullfile(fileDir,'*.ascii'));
fileNames = {fileList(:).name};

%% Get axes limits
ax = gca;
ylims = ax.YLim;
xlims = ax.XLim;

%% Add the shading to every other file
for ii=2:2:length(firsts)
    hold on
    switch timeOrSample
        case 'epsitime'
            h(ii) = fill(dataStruct.epsitime([firsts(ii),firsts(ii),lasts(ii),lasts(ii),firsts(ii)]),...
                ylims([1 2 2 1 1]),...
                [0.85 0.85 0.85],'displayname','','facealpha',0.5,'edgecolor','none');
            case 'ctdtime'            
                h(ii) = fill(dataStruct.ctdtime([firsts(ii),firsts(ii),lasts(ii),lasts(ii),firsts(ii)]),...
                ylims([1 2 2 1 1]),...
                [0.85 0.85 0.85],'displayname','','facealpha',0.5,'edgecolor','none');
        case 'epsisample'
            h(ii) = fill([firsts(ii),firsts(ii),lasts(ii),lasts(ii),firsts(ii)],...
                ylims([1 2 2 1 1]),...
                [0.85 0.85 0.85],'displayname','','facealpha',0.5,'edgecolor','none');
        case 'ctdsample'
                        h(ii) = fill([firsts(ii),firsts(ii),lasts(ii),lasts(ii),firsts(ii)],...
                ylims([1 2 2 1 1]),...
                [0.85 0.85 0.85],'displayname','','facealpha',0.5,'edgecolor','none');
        otherwise
        fprintf('Select epsitime, ctdtime, epsisample, or ctdsample for x-axis --> shade_different_files(dataStruct,timeOrSample)\n')
    end
    uistack(h(ii),'bottom')
end

%% Label every file
if strcmp(ax.YDir,'normal')
    rotAng = 90;
elseif strcmp(ax.YDir,'reverse')
    rotAng = -90;
end
switch timeOrSample
    case 'epsisample'
        for ii=1:length(firsts)
            t(ii) = text(firsts(ii)+range(xlims)/40,ylims(1)+range(ylims)/20,strrep(fileNames{ii},'_','\_'),'rotation',rotAng,'fontsize',12);
        end
    case 'ctdsample'
        for ii=1:length(firsts)
    t(ii) = text(firsts(ii)+range(xlims)/40,ylims(1)+range(ylims)/20,strrep(fileNames{ii},'_','\_'),'rotation',rotAng,'fontsize',12); 
        end
end