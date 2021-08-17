function [] = check_epsi_time(obj,saveFig)

if nargin<2
    saveFig = 0;
end

%% Get data
epsi = obj.epsi;
ctd = obj.ctd;

%% Does time ever go backwards?
    epsi_back = find(diff(epsi.epsitime)<0)+1;
    ctd_back = find(diff(ctd.ctdtime)<0)+1;
    
    %% Make samplenum arrays
    epsi.samplenum = 1:length(epsi.epsitime);
    ctd.samplenum = 1:length(ctd.ctdtime);
    
    %% Make the plot
    figure('units','inches','position',[0 0 10 13])

    % Epsitime
    ax(1) = subplot(3,1,1);
    plot(epsi.samplenum,epsi.epsitime,'.b')
    % Highlight backwards spots
    hold on
    plot(epsi.samplenum(epsi_back),epsi.epsitime(epsi_back),'dr');
    % Label new files
    shade_different_files(obj,'epsisample');
    ylabel('epsitime (sec)')
    xlabel('epsi sample num')
    
    title(strrep([obj.Meta_Data.mission ' - ' obj.Meta_Data.deployment],'_','\_'));
    
    %Ctdtime
    ax(2) = subplot(3,1,2);
    plot(ctd.samplenum,ctd.ctdtime,'.b')
    % Highlight backwards spots
    hold on
    plot(ctd.samplenum(ctd_back),ctd.ctdtime(ctd_back),'or');
    shade_different_files(obj,'ctdsample');
    ylabel('ctdtime (sec)')
    xlabel('ctd sample num')
    
    % Pressure
    ax(3) = subplot(3,1,3);
    plot(ctd.samplenum,ctd.P,'.k')
    ax(3).YDir = 'reverse';
    hold on
    shade_different_files(obj,'ctdsample');
    ylabel('pressure')
    xlabel('ctd sample num')
    
%% Save figure
if saveFig
    img = getframe(gcf);
    imwrite(img.cdata,fullfile(obj.Meta_Data.datapath,'figs/check_time.png'));
end
    
    
    

    
