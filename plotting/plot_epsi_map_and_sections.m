function [ax] = plot_epsi_map_and_sections(obj,varargin)
% [ax] = epsiPlot_map_and_sections(obj,varargin)
%
% Nicole Couto | November 2023

%% Set axes limits
if nargin<2
    limits.depth = obj.plot_properties.limits.depth;
    limits.lon = obj.plot_properties.limits.lon;
    limits.lat = obj.plot_properties.limits.lat;
    limits.temp = obj.plot_properties.limits.temp;
    limits.sal = obj.plot_properties.limits.sal;
    limits.epsilon = obj.plot_properties.limits.epsilon;
    limits.chi = obj.plot_properties.limits.chi;
else 
    limits = varargin{1};
end

%% Load the gridded profiles
data=load(fullfile(obj.Meta_Data.paths.profiles,'griddedProfiles'));
G = data.GRID;

%% Smooth bottom_depth, if you have  it
G.bottom_depth=filloutliers(G.bottom_depth,'linear');

%%  Set up figure axes
fig =  fullscreenfigure;
gap = [0.02 0.05];
marg_v = [0.08 0.05];
marg_h = [0.08  0.05];
ax(1) = subtightplot(4,3,[1,4,7,10],gap,marg_v,marg_h); %ship track
ax(2) = subtightplot(4,3,2:3,gap,marg_v,marg_h); %temperature
ax(3) = subtightplot(4,3,5:6,gap,marg_v,marg_h); %salinity
ax(4) = subtightplot(4,3,8:9,gap,marg_v,marg_h); %epsilon
ax(5) = subtightplot(4,3,11:12,gap,marg_v,marg_h); %chi
ax2pos = ax(2).Position;
drawnow
ax(6) = axes('position',[ax2pos(1),ax2pos(2)+ax2pos(4)+0.005,ax2pos(3),0.02]); %time dots

%% Plot map and ship track, colored by time
axes(ax(1))

% Contour the bathymetry
[bathy.z,bathy.lon,bathy.lat]=m_etopo2([limits.lon(1),limits.lon(2),...
                          limits.lat(1),limits.lat(2)]);
contour(bathy.lon,bathy.lat,bathy.z,'k');
hold on
contour(bathy.lon,bathy.lat,bathy.z,'k','levellist',0,'linewidth',4);

% Add the profile locations
scatter(G.longitude,G.latitude,60,G.dnum,'filled','markeredgecolor','none');
ax(1).CLim = [nanmin(G.dnum),nanmax(G.dnum)];

n_axes_map

% Add date to the title
title(['Begining ' datestr(nanmin(G.dnum),'dd-mmm-yyyy')])

%% Plot time dots above the section plots, labeled with profile number
axes(ax(6))

% Choose which profile numbers to label
if length(G.dnum)<=10
    profile_numbers = 1:length(G.dnum);
elseif length(G.dnum)>10 && length(G.dnum)<=50
    profile_numbers = 5:5:length(G.dnum);
elseif length(G.dnum)>50
    profile_numbers = 10:10:length(G.dnum);
end

scatter(G.dnum,ones(length(G.dnum),1),100,G.dnum,'filled','markeredgecolor','none');
hold on
for p=profile_numbers
    text(G.dnum(p),2.5,num2str(p),'horizontalalignment','center');
end
text(G.dnum(1),2.5,'Profile #','horizontalalignment','right');

ax(6).Visible =  'off';
ax(6).YLim = [0 2];

%% Temperature
axes(ax(2))

pcolorjw(G.dnum,G.z,G.t);
ax(2).CLim = limits.temp;
cb(2) = colorbar;
cb(2).Label.String = 'Temperature (°C)';
colormap(ax(2),cmocean('thermal'));

%% Salinity
axes(ax(3))

pcolorjw(G.dnum,G.z,G.s);
ax(3).CLim = limits.sal;
cb(3) =  colorbar;
cb(3).Label.String = 'Salinity';
colormap(ax(3),cmocean('haline'));

%% Epsilon
axes(ax(4))

pcolorjw(G.dnum,G.z,log10(G.epsilon_final));
ax(4).CLim = limits.epsilon;
cb(4) = colorbar;
cb(4).Label.String = '\epsilon (m^2 s^{-2})';
colormap(ax(4),parula);

%% Chi
axes(ax(5))
pcolorjw(G.dnum,G.z,log10(G.chi1));
ax(5).CLim = limits.chi;
cb(5) = colorbar;
cb(5).Label.String = '\chi_1 (K^2 s^{-1})';
colormap(ax(5),flipud(colorbrewer('Spectral')));

%% Add density contours
% Get density levels from the mean density profile
dens_levels_0 = n_spaced_density_levels(nanmean(G.sgth,2),G.z,20);
dens_levels_1 = rmoutliers(dens_levels_0);
dens_levels = round(dens_levels_0(1):mean(diff(dens_levels_1)):dens_levels_0(end),1);

for iAx=2:5
    axes(ax(iAx))
    hold on
    [c,ch] = contour(G.dnum,G.z,G.sgth,'k','levellist',dens_levels);
    clabel(c,ch)
end

%% Add bathymetry if you have data for at least have the profiles
if sum(~isnan(G.bottom_depth))>length(G.dnum)/2
    for iAx=2:5
        n_fill_bathy(G.dnum,G.bottom_depth);
    end
end

%% Adjust  axes properties
% Y-direction and limits
[ax(2:5).YDir] = deal('reverse');
[ax(2:5).YLim] = deal(limits.depth);

% X-limits and labels
% Adjust x-axes limits depending on how long the deployment is
if range(G.dnum)<=6/24 
    % Less than 6 hours, widen range to nearest 15 minutes
    limits.dnum = [round_time_AAA(G.dnum(1),'minute',15,'floor'),...
                   round_time_AAA(G.dnum(end),'minute',15,'ceil')];
elseif range(G.dnum)>6/24  && range(G.dnum)<1  
    % Between 6 and 24 hours, widen range to nearest 30 minutes
    limits.dnum = [round_time_AAA(G.dnum(1),'minute',30,'floor'),...
                   round_time_AAA(G.dnum(end),'minute',30,'ceil')];
elseif range(G.dnum)>1
    % Longer than 24 hours, widen range to nearest 1 hour
    limits.dnum = [round_time_AAA(G.dnum(1),'hour',1,'floor'),...
                   round_time_AAA(G.dnum(end),'hour',1,'ceil')];
end
[ax(2:5).XLim] = deal(limits.dnum);

for iAx=2:5
    datetick(ax(iAx),'x','keeplimits');  
end
[ax(2:4).XTickLabel] = deal('');

% Link x-axes
lp = linkprop([ax(2:6)],'xlim');

