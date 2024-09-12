function [ax] = fctdPlot_timeseries_spectra_tridente(obj,saveFig,ax)
% Plot timeseries of epsi, ctd, and altimeter data. For TFO2024 also
% include tridente data. I also want to plot the spectra from the last bin.
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
        too_close = obj.alt.dist<=2;
        plot(ax(8),time_array.alt(too_close),obj.alt.dst(too_close),'o','r')
    end
    
    end  
end

%% EPSI - ACC and MICROCOND PLOTS

if isclassfield(obj,'epsi') 
    if ~isempty(obj.epsi)
        
    % Acceleration z
    plot(ax(5),time_array.epsi,obj.epsi.a1_g,'.','Color',cols.a1,'LineWidth',obj.plot_properties.LineWidth,'displayname','z')

    % Acceleration x and y
    plot(ax(6),time_array.epsi,obj.epsi.a2_g,'.','Color',cols.a2,'LineWidth',obj.plot_properties.LineWidth,'displayname','x')
    hold(ax(6),'on')
    plot(ax(6),time_array.epsi,obj.epsi.a3_g,'.','Color',cols.a3,'LineWidth',obj.plot_properties.LineWidth,'displayname','y')

    % Add legends
    legend(ax(5),'location','northwest')
    legend(ax(6),'location','northwest')

    % Microconductivity
    plot(ax(4),time_array.epsi,obj.epsi.s2_volt,'.','Color',cols.ucond,'LineWidth',obj.plot_properties.LineWidth);
    
    end
end

%% FLUOROMETER

if isclassfield(obj,'fluor')

if ~isempty(obj.fluor)
    
    % Fluorometer
    plot(ax(3),time_array.fluor,obj.fluor.chla,'Color',cols.chla,'LineWidth',obj.plot_properties.LineWidth);

end

end

%% EPSI - SPECTRA PLOTS - s2, a1, a2, a3

%% GPS DATA - Print lat/lon in the corner
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

ylabel(ax(1),'T [°C]');
ylabel(ax(2),'S');
ylabel(ax(3),'chla [?]');
ylabel(ax(4),'ucond [volt]')
ylabel(ax(5),'Accel [g]');
ylabel(ax(7),'dzdt [m/s]');
ylabel(ax(8),'altimeter [m]');

[ax(:).FontSize] = deal(obj.plot_properties.FontSize);
[ax(:).FontName] = deal(obj.plot_properties.FontName);

% Right-hand y-axes
ax(6).YAxisLocation = 'right';

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
axes(ax(8))

% Set font name and size
[ax(:).FontSize] = deal(obj.plot_properties.FontSize);
[ax(:).FontName] = deal(obj.plot_properties.FontName);

% % Add tags (for tracking axes)
ax(1).Tag = 'T';
ax(2).Tag = 'S';
ax(3).Tag = 'chla';
ax(4).Tag = 'ucond';
ax(5).Tag = 'acc1';
ax(6).Tag = 'acc2_3';
ax(7).Tag = 'dzdt';
ax(8).Tag = 'alt';

if saveFig
    img = getframe(gcf);
    imwrite(img.cdata,fullfile(obj.Meta_Data.paths.figures,'epsi_ctd_alt_timeseries.png'));
    savefig(fullfile(obj.Meta_Data.paths.figures,'epsi_ctd_alt_timeseries.fig'));
end

%% SUBFUNCTIONS

function [ax] = setup_figure()

    figure('units','inches','position',[0 0 20 13])

    gap = [0.02 0.025];
    margV = [0.055 0.025];
    margH = [0.08 0.08];

    nr = 7;
    nc = 2;

    ax(1)=subtightplot(nr,nc,1,gap,margV,margH); %T
    ax(2)=subtightplot(nr,nc,3,gap,margV,margH); %S
    ax(3)=subtightplot(nr,nc,5,gap,margV,margH); %fluor
    ax(4)=subtightplot(nr,nc,7,gap,margV,margH); %microconductivity
    accAx=subtightplot(nr,nc,9,gap,margV,margH); %divide the acc plots in two
    accPos=accAx.Position;
    ax(5)=axes('position',[accPos(1) accPos(2)+accPos(4)/2 accPos(3) accPos(4)/2]); %a1
    ax(6)=axes('position',[accPos(1) accPos(2) accPos(3) accPos(4)/2]); %a2 and a3
    ax(7)=subtightplot(nr,nc,11,gap,margV,margH); %dzdt
    ax(8)=subtightplot(nr,nc,13,gap,margV,margH); %alt
    
    ax(9)=subtightplot(nr,nc,[6,8,10],[0.025 0.09],margV,[0.02 0.02]);

    delete(accAx)

function [time_array] = make_time_arrays(obj)

field_names = {'epsi','ctd','alt','vnav','gps','fluor'};

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
