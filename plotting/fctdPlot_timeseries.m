function [ax] = fctdPlot_timeseries(obj,saveFig,ax)
% Plot timeseries of epsi, ctd, and altimeter data.
%
% INPUTS
%   obj     - epsi_class object or Profile structure
%   saveFig - option to save figure
%   ax      - if running in realtime mode, input the axes handle to replace data
%
% Nicole Couto | May 2023
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
    if strcmp(ax(1).Tag,'shear')
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

%% CTD PLOTS

if isclassfield(obj,'ctd') 
    if ~isempty(obj.ctd)
        
        % Temperature
        plot(ax(1),time_array.ctd,obj.ctd.T,'.','Color',cols.T,'LineWidth',obj.plot_properties.LineWidth);
        
        % Salinity
        plot(ax(2),time_array.ctd,obj.ctd.S,'.','Color',cols.S,'LineWidth',obj.plot_properties.LineWidth);
        
        % Fall speed (dPdt)
        fall_speed = movmean(obj.ctd.dzdt,100);
        going_down = fall_speed>=0;
        coming_up = fall_speed<0;
        plot(ax(7),time_array.ctd,obj.ctd.z,'.','Color',cols.dPdt,'LineWidth',obj.plot_properties.LineWidth);
        hold(ax(7),'on')
        plot(ax(7),time_array.ctd(going_down),obj.ctd.z(going_down),'.','Color',[0.2157 0.4941 0.7216],'LineWidth',obj.plot_properties.LineWidth);
        plot(ax(7),time_array.ctd(coming_up),obj.ctd.z(coming_up),'.','Color',[0.9686 0.5059 0.7490],'LineWidth',obj.plot_properties.LineWidth);
        
        % Add arrows showing up and down
        annotation(gcf,'arrow',[0.94 0.94], [0.42 0.38],'Color',[0.2157 0.4941 0.7216]);
        annotation(gcf,'arrow',[0.96 0.96], [0.38 0.42],'Color',[0.9686 0.5059 0.7490]);

        % Alitmeter
        if isclassfield(obj,'alt') && ~isempty(obj.alt)
            plot(ax(8),time_array.alt,obj.alt.dst,'.','Color',cols.alt)
            hold(ax(8),'on');
            too_close = obj.alt.dst<=2;
            if sum(too_close)>0
                plot(ax(8),time_array.alt(too_close),obj.alt.dst(too_close),'o','r');
            end
        end
    
    end  
end

%% VNAV PLOTS

if isclassfield(obj,'vnav') 
    if ~isempty(obj.vnav)
        
    % Gyro z
    plot(ax(5),time_array.vnav,obj.vnav.gyro(:,3),'.','Color',cols.gyro3,'LineWidth',obj.plot_properties.LineWidth,'displayname','z')

    % Gyro x
    plot(ax(6),time_array.vnav,obj.vnav.gyro(:,1),'.','Color',cols.gyro1,'LineWidth',obj.plot_properties.LineWidth,'displayname','x')
    hold(ax(6),'on')
    plot(ax(6),time_array.vnav,obj.vnav.gyro(:,2),'.','Color',cols.gyro2,'LineWidth',obj.plot_properties.LineWidth,'displayname','y')

    % Add legends
    legend(ax(5),'location','northwest')
    legend(ax(6),'location','northwest')
    
    end
end

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
        datetick(ax(10),'x','HH:MM:SS','keepticks')
        ax(10).XLabel.String = 'HH:MM:SS';
    catch
    end
elseif ~replaceData
    [ax(:).XLim] = deal([nanmin(obj.epsi.dnum),nanmax(obj.epsi.dnum)]);
    try
        datetick(ax(10),'x','MM:SS','keepticks')
        ax(10).XLabel.String = 'MM:SS';
    catch
    end
else
    ax(10).XLabel.String = 'epsitime (s)';
end

% Title, ylabels, and font size
if isfield(obj,'Meta_Data')
    title(ax(1),sprintf('%s-%s-%s',strrep(obj.Meta_Data.mission,'_','\_'),...
        strrep(obj.Meta_Data.vehicle_name,'_','\_'),...
        strrep(obj.Meta_Data.deployment,'_','\_')));
end

% Labels
[ax(1:9).XTickLabel]=deal('');

ylabel(ax(1),'T [°C]');
ylabel(ax(2),'S');
ylabel(ax(3),'Accel [g]');
ylabel(ax(5),'Gyro []');
ylabel(ax(7),'z [m]');
ylabel(ax(9),'altimeter');

[ax(:).FontSize] = deal(obj.plot_properties.FontSize);
[ax(:).FontName] = deal(obj.plot_properties.FontName);

% Right-hand y-axes
ax(4).YAxisLocation = 'right';
ax(6).YAxisLocation = 'right';
ax(10).YAxisLocation = 'right';
ax(10).Color = 'none';
ax(10).YLim = [0 30];
ax(10).YGrid = 'on';

% Flip pressure axis
ax(8).YDir = 'reverse';

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

% Bring alt axes to the front
axes(ax(10))

% Set font name and size
[ax(:).FontSize] = deal(obj.plot_properties.FontSize);
[ax(:).FontName] = deal(obj.plot_properties.FontName);

% % Add tags (for tracking axes)
ax(1).Tag = 'shear';
ax(2).Tag = 'fpo7';
ax(3).Tag = 'a1';
ax(4).Tag = 'a2_a3';
ax(5).Tag = 'gyro3';
ax(6).Tag = 'gyro1_2';
ax(7).Tag = 'compass';
ax(8).Tag = 'dPdt';
ax(9).Tag = 'P';
ax(10).Tag = 'alt';

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
    gyroAx=subtightplot(7,1,4,gap,margV,margH); %divide the acc plots in two
    gyroPos=gyroAx.Position;
    ax(5)=axes('position',[gyroPos(1) gyroPos(2)+gyroPos(4)/2 gyroPos(3) gyroPos(4)/2]); %gyro z
    ax(6)=axes('position',[gyroPos(1) gyroPos(2) gyroPos(3) gyroPos(4)/2]); %gyro x and y
    ax(7)=subtightplot(7,1,5,gap,margV,margH); %compass
    ax(8)=subtightplot(7,1,6,gap,margV,margH); %dPdt
    ax(9)=subtightplot(7,1,7,gap,margV,margH); %P
    ax(10)=axes('position',ax(9).Position);     %alt
    
    delete(accAx)
    delete(gyroAx)
    
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
