function [ax] = epsiPlot_timeseries(obj,saveFig,ax)
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
nSec = 1*60;
nDay = nSec/(3600*24);

%% Set plot properties if you don't have them
if ~isclassfield(obj,'plot_properties')
    obj.plot_properties = epsiSetup_set_plot_properties;
end

%% Define time axes for the plot
time_array = make_time_arrays(obj);

%% Set up axes
if ~replaceData
    ax = setup_figure;
elseif replaceData
    if strcmp(ax(1).Tag,'fp07')
        for iAx=1:length(ax)
            ax(iAx).NextPlot = 'replace';
        end
        % Delete the lat/lon textbox, if it exists
        fig = ax(1).Parent;
        tb = findall(fig,'type','textboxshape');
        delete(tb)
        
        % Delete any other figures (like the blank one set up when you call
        % this in a timer function)
        all_figs = findobj('type','figure');
        for f=1:length(all_figs)
            fig_axes = findobj(all_figs(f),'type','axes');
            old = 0;
            for a=1:length(fig_axes)
                if isempty(fig_axes(a).Tag) %I give all the axes a tag name, so if this is empty, it's an old empty figure
                    old = 1;
                end
            end
            if old==1
                delete(all_figs(f))
            end
        end
        
    else
        ax = setup_figure;
    end 
end
cols = obj.plot_properties.Colors;

%% EPSI PLOTS

if isclassfield(obj,'epsi') && ~isempty(obj.epsi)
    
    % t1 and t2 - differences from mean
    plot(ax(1),time_array.epsi,obj.epsi.t1_volt-nanmean(obj.epsi.t1_volt),'.','Color',cols.t1,'LineWidth',obj.plot_properties.LineWidth,'displayname',sprintf('t1 - %1.1f',nanmean(obj.epsi.t1_volt)));
    hold(ax(1),'on')
    plot(ax(1),time_array.epsi,obj.epsi.t2_volt-nanmean(obj.epsi.t2_volt),'.','Color',cols.t2,'LineWidth',obj.plot_properties.LineWidth,'displayname',sprintf('t2 - %1.1f',nanmean(obj.epsi.t2_volt)));
    hold(ax(1),'off')
    
%     % t1 and t2
%     plot(ax(1),time_array.epsi,obj.epsi.t1_volt,'.','Color',cols.t1,'LineWidth',obj.plot_properties.LineWidth,'displayname','t1');
%     hold(ax(1),'on')
%     plot(ax(1),time_array.epsi,obj.epsi.t2_volt,'.','Color',cols.t2,'LineWidth',obj.plot_properties.LineWidth,'displayname','t2');
%    
    
    % s1 and s2
    plot(ax(2),time_array.epsi,obj.epsi.s1_volt,'.','Color',cols.s1,'LineWidth',obj.plot_properties.LineWidth,'displayname',sprintf('s1 -rms %1.1f',rms(obj.epsi.s1_volt)))
    hold(ax(2),'on')
    plot(ax(2),time_array.epsi,obj.epsi.s2_volt,'.','Color',cols.s2,'LineWidth',obj.plot_properties.LineWidth,'displayname',sprintf('s2 - rms%1.1f',rms(obj.epsi.s2_volt)))

    % a1
    plot(ax(3),time_array.epsi,obj.epsi.a1_g,'.','Color',cols.a1,'LineWidth',obj.plot_properties.LineWidth,'displayname','a1')
    
    % a2 and a3
    plot(ax(4),time_array.epsi,obj.epsi.a2_g,'.','Color',cols.a2,'LineWidth',obj.plot_properties.LineWidth,'displayname','a2')
    hold(ax(4),'on')
    plot(ax(4),time_array.epsi,obj.epsi.a3_g,'.','Color',cols.a3,'LineWidth',obj.plot_properties.LineWidth,'displayname','a3')

    % Add legends
    legend(ax(1),'location','northwest')
    legend(ax(2),'location','northwest')
    legend(ax(3),'location','northwest')
    legend(ax(4),'location','northwest')
end

%% CTD PLOTS

if isclassfield(obj,'ctd') 
    if ~isempty(obj.ctd)

        % Temperature
        plot(ax(5),time_array.ctd,obj.ctd.T,'.','Color',cols.T,'LineWidth',obj.plot_properties.LineWidth);
        
        % Salinity
        plot(ax(6),time_array.ctd,obj.ctd.S,'.','Color',cols.S,'LineWidth',obj.plot_properties.LineWidth);

        % Depth and fall speed
       fall_speed = movmean(obj.ctd.dzdt,100,'omitmissing');
        going_down = fall_speed>=0;
        coming_up = fall_speed<0;
        plot(ax(7),time_array.ctd,fall_speed,'.','Color',cols.dPdt,'LineWidth',obj.plot_properties.LineWidth);
        hold(ax(7),'on')
        plot(ax(7),time_array.ctd(going_down),fall_speed(going_down),'.','Color',[0.2157 0.4941 0.7216],'LineWidth',obj.plot_properties.LineWidth);
        plot(ax(7),time_array.ctd(coming_up),fall_speed(coming_up),'.','Color',[0.9686 0.5059 0.7490],'LineWidth',obj.plot_properties.LineWidth);
        
        % Add arrows showing up and down
        annotation(gcf,'arrow',[0.94 0.94], [0.29 0.25],'Color',[0.2157 0.4941 0.7216]);
        annotation(gcf,'arrow',[0.96 0.96], [0.25 0.29],'Color',[0.9686 0.5059 0.7490]);

        % Alitmeter
        if isclassfield(obj,'alt') && ~isempty(obj.alt)
            plot(ax(8),time_array.alt,obj.alt.dst,'.','Color',cols.alt)
            hold(ax(8),'on');
            too_close = obj.alt.dst<=2;
            if sum(too_close)>0
                plot(ax(8),time_array.alt(too_close),obj.alt.dst(too_close),'or');
            end
        end
    
    end  
end

% %% VNAV PLOTS
% 
% if isclassfield(obj,'vnav') 
%     if ~isempty(obj.vnav)
% 
%     % Gyro z
%     plot(ax(5),time_array.vnav,obj.vnav.gyro(:,3),'.','Color',cols.gyro3,'LineWidth',obj.plot_properties.LineWidth,'displayname','z')
% 
%     % Gyro x
%     plot(ax(6),time_array.vnav,obj.vnav.gyro(:,1),'.','Color',cols.gyro1,'LineWidth',obj.plot_properties.LineWidth,'displayname','x')
%     hold(ax(6),'on')
%     plot(ax(6),time_array.vnav,obj.vnav.gyro(:,2),'.','Color',cols.gyro2,'LineWidth',obj.plot_properties.LineWidth,'displayname','y')
% 
%     % Compass
%     plot(ax(7),time_array.vnav,obj.vnav.compass(:,1),'.','Color',cols.compass1,'LineWidth',obj.plot_properties.LineWidth);
% 
%     % Add legends
%     legend(ax(5),'location','northwest')
%     legend(ax(6),'location','northwest')
% 
%     end
% end

%% GPS DATA
if isclassfield(obj,'gps') 
    if ~isempty(obj.gps)
        % Find the last non-nan lat/lon and plot a few decimals
        not_nan_lat = obj.gps.latitude(~isnan(obj.gps.latitude) & ~isnan(obj.gps.longitude));
        not_nan_lon = obj.gps.longitude(~isnan(obj.gps.latitude) & ~isnan(obj.gps.longitude));
        if ~isempty(not_nan_lat)
            lat = not_nan_lat(end);
            lon = not_nan_lon(end);
        else
            lat = nan;
            lon = nan;
        end
        
        annotation(gcf,'textbox',...
            [0.8093 0.9754 0.2212 0.0288],...
            'String',{sprintf('lat = %03.4f',lat),sprintf('lon = %03.4f',lon)},...
            'LineStyle','none',...
            'FontSize',14,...
            'FontName',obj.plot_properties.FontName,...
            'FitBoxToText','off');
    end
end


%% AXES

% Time axes label and limits with 10-sec tick marks
if isfield(obj.epsi,'dnum') && ~all(isnan(obj.epsi.dnum)) && replaceData
    sec10 = 10/(3600*24);
    % If plotting in realtime, limit view
    [ax(:).XTick] = deal(fliplr(nanmax(obj.epsi.dnum):-sec10:nanmax(obj.epsi.dnum)-nDay));
    [ax(:).XLim] = deal([nanmax(obj.epsi.dnum)-nDay,nanmax(obj.epsi.dnum)]);
    try
        datetick(ax(8),'x','HH:MM:SS','keepticks')
        ax(8).XLabel.String = 'HH:MM:SS';
    catch
    end
elseif ~replaceData
    [ax(:).XLim] = deal([nanmin(obj.epsi.dnum),nanmax(obj.epsi.dnum)]);
    try
        datetick(ax(8),'x','MM:SS','keepticks')
        ax(8).XLabel.String = 'MM:SS';
    catch
    end
else
    ax(8).XLabel.String = 'epsitime (s)';
end

% Title, ylabels, and font size
if isfield(obj,'Meta_Data')
    title(ax(1),sprintf('%s-%s-%s',strrep(obj.Meta_Data.mission,'_','\_'),...
        strrep(obj.Meta_Data.vehicle_name,'_','\_'),...
        strrep(obj.Meta_Data.deployment,'_','\_')));
end

% Labels
[ax(1:7).XTickLabel]=deal('');

ylabel(ax(1),'FPO7 [Volt]');
ylabel(ax(2),'Shear [Volt]');
ylabel(ax(3),'Accel [g]');
ylabel(ax(5),'T [°C]');
ylabel(ax(6),'S');
ylabel(ax(7),'dz/dt [m/s]'); % ALB TODO add z also 
ylabel(ax(8),'altimeter');

[ax(:).FontSize] = deal(obj.plot_properties.FontSize);
[ax(:).FontName] = deal(obj.plot_properties.FontName);

% Right-hand y-axes
ax(4).YAxisLocation = 'right';
ax(8).Color = 'none';
ax(8).YLim = [0 30];
ax(8).YGrid = 'on';

% Flip pressure axis
ax(7).YDir = 'reverse';

% Link x-axes and add grid
linkprop([ax(:)],'xlim');
[ax(:).XGrid] = deal('on');

% Y-Limits
% for iAx=1:length(ax)
%     dataMed = nanmedian([ax(iAx).Children(:).YData]);
%     dataStd = nanstd([ax(iAx).Children(:).YData]);
%     try %Try to reset the y-limits. Sometimes it fails if all data is nan
%     ax(iAx).YLim = [dataMed-1.5*dataStd,dataMed+1.5*dataStd];
%     catch
%     end
% end

% Set font name and size
[ax(:).FontSize] = deal(obj.plot_properties.FontSize);
[ax(:).FontName] = deal(obj.plot_properties.FontName);

% % Add tags (for tracking axes)
ax(1).Tag = 'fp07';
ax(2).Tag = 'shear';
ax(3).Tag = 'a1';
ax(4).Tag = 'a2_a3';
ax(5).Tag = 'T';
ax(6).Tag = 'S';
ax(7).Tag = 'P and dzdt';
ax(8).Tag = 'alt';

if saveFig
    img = getframe(gcf);
    imwrite(img.cdata,fullfile(obj.Meta_Data.paths.figures,'epsi_ctd_alt_timeseries.png'));
    savefig(fullfile(obj.Meta_Data.paths.figures,'epsi_ctd_alt_timeseries.fig'));
end

%% SUBFUNCTIONS

function [ax] = setup_figure()

figure('units','inches','position',[0 0 10 13])


    gap = [0.02 0.025];
    margV = [0.045 0.035];
    margH = [0.12 0.08];

    ax(1)=subtightplot(7,1,1,gap,margV,margH); %shear
    ax(2)=subtightplot(7,1,2,gap,margV,margH); %fpo7
    accAx=subtightplot(7,1,3,gap,margV,margH); %divide the acc plots in two
    accPos=accAx.Position;
    ax(3)=axes('position',[accPos(1) accPos(2)+accPos(4)/2 accPos(3) accPos(4)/2]); %a1
    ax(4)=axes('position',[accPos(1) accPos(2) accPos(3) accPos(4)/2]); %a2 and a3
    ax(5)=subtightplot(7,1,4,gap,margV,margH); %T
    ax(6)=subtightplot(7,1,5,gap,margV,margH); %S
    ax(7)=subtightplot(7,1,6,gap,margV,margH); %P and dzdt
    ax(8)=subtightplot(7,1,7,gap,margV,margH); %alt
    
    delete(accAx)
    %delete(gyroAx)
    
%     gap = [0.025 0.025];
%     margV = [0.08 0.05];
%     margH = [0.1 0.1];
% 
%     ax(1)=subtightplot(6,1,1,gap,margV,margH); %shear
%     ax(2)=subtightplot(6,1,2,gap,margV,margH); %fpo7
%     accAx=subtightplot(6,1,3,gap,margV,margH); %divide the acc plots in two
%     accPos=accAx.Position;
%     ax(3)=axes('position',[accPos(1) accPos(2)+accPos(4)/2 accPos(3) accPos(4)/2]); %a1
%     ax(4)=axes('position',[accPos(1) accPos(2) accPos(3) accPos(4)/2]); %a2 and a3
%     delete(accAx)
%     ax(5)=subtightplot(6,1,4,gap,margV,margH); %P
%     ax(6)=axes('position',ax(5).Position);     %alt
%     ax(7)=subtightplot(6,1,5,gap,margV,margH); %dPdt
%     ax(8)=subtightplot(6,1,6,gap,margV,margH); %T
%     ax(9)=axes('position',ax(8).Position);     %S



function [time_array] = make_time_arrays(obj)

field_names = {'epsi','ctd','alt','vnav','gps'};

for iField = 1:length(field_names)

classCondition = isclassfield(obj,field_names{iField}) && ~isempty(obj.(field_names{iField}));
if classCondition
    structCondition = ~isempty(obj.(field_names{iField}));
    if structCondition
        % Do you have dnum or seconds?
        if isfield(obj.(field_names{iField}),'dnum')
            time_array.(field_names{iField}) = obj.(field_names{iField}).dnum;
        else
            time_array.(field_names{iField}) = obj.(field_names{iField}).time_s;
        end
    end
end

end
