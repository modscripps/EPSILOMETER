function [] = shade_different_files(epsi,timeOrSamplenum)

% Add shaded boxes to a plot to denote different files
%
% USAGE:
%   [] = shade_different_files(epsi,timeOrSamplenum)
%
% INPUTS:
%   epsi            = epsi structure
%   timeOrSamplenum = 'time' if x-axis is time, or 'sample' if x-axis is
%                      sample number


%% Find where filenum changes
dF = find(diff(epsi.fileNum)~=0);
firsts = [1;dF];
lasts = [dF+1;length(epsi.fileNum)];

%% Add the shading to every other file
for ii=2:2:length(firsts)
    ax = gca;
    hold on
    ylims = ax.YLim;
    switch timeOrSamplenum
        case 'time'
            h(ii) = fill(epsi.epsitime([firsts(ii),firsts(ii),lasts(ii),lasts(ii),firsts(ii)]),...
                ylims([1 2 2 1 1]),...
                [0.85 0.85 0.85],'displayname','','facealpha',0.5,'edgecolor','none');
        case 'samplenum'
            h(ii) = fill([firsts(ii),firsts(ii),lasts(ii),lasts(ii),firsts(ii)],...
                ylims([1 2 2 1 1]),...
                [0.85 0.85 0.85],'displayname','','facealpha',0.5,'edgecolor','none');
        otherwise
        fprintf('Choose time or samplenum for x-axis\n')
    end
    uistack(h(ii),'bottom')
end