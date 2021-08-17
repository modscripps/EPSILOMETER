function [ax] = epsiPlot_epsi_ctd_alt_timeseries(obj,saveFig,ax)
% Plot timeseries of epsi, ctd, and altimeter data.
%
% INPUTS
%   obj     - epsi_class object or Profile structure
%   saveFig - option to save figure
%   ax      - if running in realtime mode, input the axes handle to replace data
%
% Nicole Couto | June 2021
% -------------------------------------------------------------------------
if nargin==3
    replaceData=1;
elseif nargin<3
    replaceData=0;
    if nargin<2
        saveFig=0;
    end
end

%% If you're running in realtime, view the most recent nSec
nSec = 2*60;
nDay = nSec/(3600*24);

%% Set plot properties if you don't have them
if ~isfield(obj,'plot_properties')
    obj.plot_properties = set_epsi_plot_properties;
end

%% Define time axes for the plot
% Epsi
classCondition = isclassfield(obj,'epsi') && ~isempty(obj.epsi);
structCondition = ~isempty(obj.epsi);
if classCondition || structCondition
    % Do you have dnum or seconds?
    if isfield(obj.epsi,'dnum')
        timeEpsi = obj.epsi.dnum;
        nTime = nDay;
    else
        timeEpsi = obj.epsi.time_s;
        nTime = nSec;
    end
end

% Ctd
classCondition = isclassfield(obj,'ctd') && ~isempty(obj.ctd);
structCondition = ~isempty(obj.ctd);
if classCondition || structCondition
     % Do you have dnum or seconds?
    if isfield(obj.ctd,'dnum')
        timeCtd = obj.ctd.dnum;
        nTime = nDay;
    else
        timeCtd = obj.ctd.time_s;
        nTime = nSec;  
    end
end

% Altimeter
classCondition = isclassfield(obj,'alt') && ~isempty(obj.alt);
structCondition = ~isempty(obj.alt);
if classCondition || structCondition
    % Do you have dnum or seconds?
    if isfield(obj.alt,'dnum')
        timeAlt = obj.alt.dnum;
        nTime = nDay;
    else
        timeAlt = obj.alt.time_s;
        nTime = nSec;
    end
end

%% Set up axes
if ~replaceData
    figure('units','inches','position',[0 0 10 13])
    
    gap = [0.025 0.025];
    margV = [0.08 0.05];
    margH = [0.1 0.1];
    
    ax(1)=subtightplot(6,1,1,gap,margV,margH); %shear
    ax(2)=subtightplot(6,1,2,gap,margV,margH); %fpo7
    accAx=subtightplot(6,1,3,gap,margV,margH); %divide the acc plots in two
    accPos=accAx.Position;
    ax(3)=axes('position',[accPos(1) accPos(2)+accPos(4)/2 accPos(3) accPos(4)/2]); %a1
    ax(4)=axes('position',[accPos(1) accPos(2) accPos(3) accPos(4)/2]); %a2 and a3
    delete(accAx)
    ax(5)=subtightplot(6,1,4,gap,margV,margH); %P
    ax(6)=axes('position',ax(5).Position);     %alt
    ax(7)=subtightplot(6,1,5,gap,margV,margH); %dPdt
    ax(8)=subtightplot(6,1,6,gap,margV,margH); %T
    ax(9)=axes('position',ax(8).Position);     %S
    
elseif replaceData
    if strcmp(ax(1).Tag,'shear')
        for iAx=1:length(ax)
            ax(iAx).NextPlot = 'replace';
        end
    else
        figure('units','inches','position',[0 0 10 13])
        
        gap = [0.025 0.025];
        margV = [0.08 0.05];
        margH = [0.1 0.1];
        
        ax(1)=subtightplot(6,1,1,gap,margV,margH); %shear
        ax(2)=subtightplot(6,1,2,gap,margV,margH); %fpo7
        accAx=subtightplot(6,1,3,gap,margV,margH); %divide the acc plots in two
        accPos=accAx.Position;
        ax(3)=axes('position',[accPos(1) accPos(2)+accPos(4)/2 accPos(3) accPos(4)/2]); %a1
        ax(4)=axes('position',[accPos(1) accPos(2) accPos(3) accPos(4)/2]); %a2 and a3
        delete(accAx)
        ax(5)=subtightplot(6,1,4,gap,margV,margH); %P
        ax(6)=axes('position',ax(5).Position);     %alt
        ax(7)=subtightplot(6,1,5,gap,margV,margH); %dPdt
        ax(8)=subtightplot(6,1,6,gap,margV,margH); %T
        ax(9)=axes('position',ax(8).Position);     %S
        
    end
end
cols = obj.plot_properties.Colors;

%% EPSI PLOTS

if isclassfield(obj,'epsi') && ~isempty(obj.epsi)
     
    
%     plot(ax(1),timeEpsi,obj.epsi.t1_volt - nanmean(obj.epsi.t1_volt),'.','Color',cols.t1,'LineWidth',obj.plot_properties.LineWidth,'displayname','t1 (mean removed)')
%     hold(ax(1),'on')
%     plot(ax(1),timeEpsi,obj.epsi.t2_volt - nanmean(obj.epsi.t2_volt),'.',...'Color',cols.t2,'LineWidth',obj.plot_properties.LineWidth,'displayname','t2 (mean removed)')
    
    plot(ax(1),timeEpsi,obj.epsi.t1_volt,'.','Color',cols.t1,'LineWidth',obj.plot_properties.LineWidth,'displayname','t1')
    hold(ax(1),'on')
    plot(ax(1),timeEpsi,obj.epsi.t2_volt,'.','Color',cols.t2,'LineWidth',obj.plot_properties.LineWidth,'displayname','t2')
    
    plot(ax(2),timeEpsi,obj.epsi.s1_volt,'.','Color',cols.s1,'LineWidth',obj.plot_properties.LineWidth,'displayname','s1')
    hold(ax(2),'on')
    plot(ax(2),timeEpsi,obj.epsi.s2_volt,'.','Color',cols.s2,'LineWidth',obj.plot_properties.LineWidth,'displayname','s2')
    
    plot(ax(3),timeEpsi,obj.epsi.a1_g,'.','Color',cols.a1,'LineWidth',obj.plot_properties.LineWidth,'displayname','a1')
    
    plot(ax(4),timeEpsi,obj.epsi.a2_g,'.','Color',cols.a2,'LineWidth',obj.plot_properties.LineWidth,'displayname','a2')
    hold(ax(4),'on')
    plot(ax(4),timeEpsi,obj.epsi.a3_g,'.','Color',cols.a3,'LineWidth',obj.plot_properties.LineWidth,'displayname','a3')
    
    legend(ax(1),'location','northwest')
    legend(ax(2),'location','northwest')
    legend(ax(3),'location','northwest')
    legend(ax(4),'location','northwest')
    
    % Title, ylabels, and font size
    if isfield(obj,'Meta_Data')
    title(ax(1),sprintf('%s-%s-%s',strrep(obj.Meta_Data.mission,'_','\_'),...
        strrep(obj.Meta_Data.vehicle_name,'_','\_'),...
        strrep(obj.Meta_Data.deployment,'_','\_')));
    end
    ylabel(ax(1),'FPO7 [Volt]','Fontsize',obj.plot_properties.FontSize,'FontName',obj.plot_properties.FontName)
    ylabel(ax(2),'Shear [Volt]','Fontsize',obj.plot_properties.FontSize,'FontName',obj.plot_properties.FontName)
    ylabel(ax(3),'Accel [g]','Fontsize',obj.plot_properties.FontSize,'FontName',obj.plot_properties.FontName)
    %ax(3).XLabel.String = 'seconds';
    set(ax(1),'Fontsize',obj.plot_properties.FontSize,'FontName',obj.plot_properties.FontName)
    set(ax(2),'Fontsize',obj.plot_properties.FontSize,'FontName',obj.plot_properties.FontName)
    set(ax(3),'Fontsize',obj.plot_properties.FontSize,'FontName',obj.plot_properties.FontName)
    
end

%% CTD PLOTS

if isclassfield(obj,'ctd') && ~isempty(obj.ctd)
    
    % Pressure and altimiter
    plot(ax(5),timeCtd,obj.ctd.P,'.','Color',cols.P)
    if isclassfield(obj,'alt') && ~isempty(obj.alt)
        plot(ax(6),timeAlt,obj.alt.dst,'.','Color',cols.alt)
    end
    
    % Fall speed (dPdt)
    plot(ax(7),timeCtd,movmean(obj.ctd.dPdt,100),'.','Color',cols.dPdt,'LineWidth',obj.plot_properties.LineWidth)
    
    % Temperature and salinity
    plot(ax(8),timeCtd,obj.ctd.T,'.','Color',cols.T,'LineWidth',obj.plot_properties.LineWidth)
    plot(ax(9),timeCtd,real(obj.ctd.S),'.','Color',cols.S,'LineWidth',obj.plot_properties.LineWidth)
    
    % Time axes label and limits with 30-sec tick marks
    if isfield(obj.ctd,'dnum') && ~all(isnan(obj.ctd.dnum))
        sec30 = 30/(3600*24);
            % Find the closest 30 second or 60 second interval
            %xMaxIdx = max(find(second(obj.ctd.dnum)==30 | second(obj.ctd.dnum)==0));
            %             if isempty(xMaxIdx)
            %                 [ax(:).XTick] = deal(fliplr(obj.ctd.dnum(xMaxIdx):-sec30:nanmax(obj.ctd.dnum)-nDay));
            %             else
            % If plotting in realtime, limit view
                [ax(:).XTick] = deal(fliplr(nanmax(obj.ctd.dnum):-sec30:nanmax(obj.ctd.dnum)-nDay));
                %end
                [ax(:).XLim] = deal([nanmax(obj.ctd.dnum)-nDay,nanmax(obj.ctd.dnum)]);
            %datetick(ax(8),'x','HH:MM:SS','keeplimits','keepticks') %Hidden behind ax(9) but we need the same ticks
            try
                datetick(ax(9),'x','HH:MM:SS','keepticks')
            catch
            end
    else
        ax(9).XLabel.String = 'ctdtime (s)';
    end
    ax(9).XTickLabelRotation = 45;
end

%% AXES

% Right-hand y-axes
ax(6).YAxisLocation = 'right';
ax(9).YAxisLocation = 'right';
ax(6).Color = 'none';
ax(9).Color = 'none';

% Axes colors
ax(5).YColor = cols.P;
ax(6).YColor = cols.alt;
ax(8).YColor = cols.T;
ax(9).YColor = cols.S;

% Flip pressure axis
ax(5).YDir = 'reverse';

    % Link x-axes and add grid
linkprop([ax(:)],'xlim');
[ax(:).XGrid] = deal('on');

% Labels
[ax(1:8).XTickLabel]=deal('');
ylabel(ax(5),'P [dbar]','Fontsize',obj.plot_properties.FontSize,'FontName',obj.plot_properties.FontName)
ylabel(ax(6),'altimeter','Fontsize',obj.plot_properties.FontSize,'FontName',obj.plot_properties.FontName)
ylabel(ax(7),'dPdt [dbar/s]','Fontsize',obj.plot_properties.FontSize,'FontName',obj.plot_properties.FontName)
ylabel(ax(8),'T [°C]','Fontsize',obj.plot_properties.FontSize,'FontName',obj.plot_properties.FontName)
ylabel(ax(9),'S [psu]','Fontsize',obj.plot_properties.FontSize,'FontName',obj.plot_properties.FontName)

% Y-Limits
% for iAx=1:length(ax)
%     dataMed = nanmedian([ax(iAx).Children(:).YData]);
%     dataStd = nanstd([ax(iAx).Children(:).YData]);
%     try %Try to reset the y-limits. Sometimes it fails if all data is nan
%     ax(iAx).YLim = [dataMed-1.5*dataStd,dataMed+1.5*dataStd];
%     catch
%     end
% end

% Bring alt and salinity axes to the front
axes(ax(6))
axes(ax(9))

[ax(:).FontSize] = deal(obj.plot_properties.FontSize);
[ax(:).FontName] = deal(obj.plot_properties.FontName);

% % Add tags (for tracking axes)
ax(1).Tag = 'shear';
ax(2).Tag = 'fpo7';
ax(3).Tag = 'a1';
ax(4).Tag = 'a2_a3';
ax(5).Tag = 'P';
ax(6).Tag = 'alt';
ax(7).Tag = 'dPdt';
ax(8).Tag = 'T';
ax(9).Tag = 'S';

if saveFig
    img = getframe(gcf);
    imwrite(img.cdata,fullfile(obj.Meta_Data.datapath,'figs/epsi_ctd_alt_timeseries.png'));
    savefig(fullfile(obj.Meta_Data.datapath,'figs/epsi_ctd_timeseries_alt.fig'));
end
