% plot_fctd_sections
fctd_mat_dir = fullfile(ec.Meta_Data.paths.data,'fctd_mat');

%% First, concatenate the individual fctd files in the deployment directory
[FCTDall,FCTDgrid] = concatenate_and_grid_fctd(fctd_mat_dir);
%% Plot some stuff
clf
%%ax = [];
fig = figure(1);
fig.Units = 'normalized';
fig.Position = [0 0 0.6 0.9];
zlim = [0 80];
clim_temp = [25.6 30.2]; 
clim_sal = [32.6 33.9];
 clim_chi = [0.2 1];
levels_dens = [19:1:25.5 26:0.2:27.7 27.71:0.01:27.8];

% which data to plot? How about the most recent 1 day
iplot=find(FCTDgrid.time>FCTDgrid.time-1); 
iplot=iplot(1:2:end); % just the down-casts till we correct the hysteresis later


% Temperature
ax(1) = subtightplot(4,1,1);
pcolorjw(FCTDgrid.time(iplot),FCTDgrid.depth,FCTDgrid.temperature(:,iplot));
ax(1).CLim = [clim_temp(1) clim_temp(2)];
cb(1) = colorbar;
colormap(ax(1),lansey)
cb(1).Label.String = 'Temperature';

% Salinity
ax(2) = subtightplot(4,1,2);
pcolorjw(FCTDgrid.time(iplot),FCTDgrid.depth,real(FCTDgrid.salinity(:,iplot)));
ax(2).CLim = [clim_sal(1) clim_sal(2)];
cb(2) = colorbar;
colormap(ax(2),cmocean('delta'))
cb(2).Label.String = 'Salinity';
set(cb(2),'ydir','reverse');

% Chlorophyll
ax(3) = subtightplot(4,1,3);
pcolorjw(FCTDgrid.time(iplot),FCTDgrid.depth,((FCTDgrid.chla(:,iplot))/2^16-0.5)*500.0);
%pcolorjw(FCTDgrid.time,FCTDgrid.depth,((FCTDgrid.bb)/2^16-0.5)*10.0);
%pcolorjw(FCTDgrid.time,FCTDgrid.depth,((FCTDgrid.fDOM)/2^16-0.5)*1000.0);
 ax(3).CLim = [clim_chi(1) clim_chi(2)];
cb(3) = colorbar;
cb(3).Label.String = 'Chla';


% N^2
ax(4) = subtightplot(4,1,4);
%pcolorjw(FCTDgrid.time,FCTDgrid.depth,((FCTDgrid.chla)/2^16-0.5)*500.0);
[bfrq,vort,p_ave] = sw_bfrq(FCTDgrid.salinity,FCTDgrid.temperature,FCTDgrid.pressure,mean(FCTDgrid.latitude,'omitmissing'));
pcolorjw(FCTDgrid.time(iplot),p_ave(:,1),real(log10(bfrq(:,iplot))));
colormap(ax(4),cmocean('speed'))
ax(4).CLim = ([-5.5 -2.5]);
cb(4) = colorbar;
cb(4).Label.String = 'N^2';

% Add density contours and datetick
for iAx=1:4
   axes(ax(iAx))
   hold(ax(iAx),'on')
   [c,ch] = contour(FCTDgrid.time,FCTDgrid.depth,real(FCTDgrid.density-1000),['w'],'levellist',levels_dens);
 %contour(FCTDgrid.time,FCTDgrid.depth,FCTDgrid.temperature,'m','levellist',13);
 %  contour(FCTDgrid.time,FCTDgrid.depth,FCTDgrid.temperature,'c','levellist',15);
   clabel(c,ch);
   
   datetick(ax(iAx),'x','HH:MM','keeplimits')
end

% Depth axes
[ax(1:4).YLim] = deal([zlim(1) zlim(2)]);
[ax(1:4).YDir] = deal('reverse');

%MHA hack
%[ax(1:2).YLim] = deal([zlim(1) 500]);



