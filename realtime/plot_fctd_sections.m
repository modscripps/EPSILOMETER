% plot_fctd_sections
fctd_mat_dir = fullfile(ec.Meta_Data.paths.data,'fctd_mat');

%% First, concatenate the individual fctd files in the deployment directory
[FCTDall,FCTDgrid] = concatenate_and_grid_fctd(fctd_mat_dir);

%% Plot some stuff
clf

fig = figure(1);
fig.Position = [0 0 960 984];
zlim = [0 1000];
clim_temp = [8 26];
clim_sal = [34.5 37];
clim_chi = [-10 -8];
levels_dens = [1025 1026:0.1:1027.7 1027.71:0.01:1027.8];

% Temperature
ax(1) = subtightplot(3,1,1);
pcolorjw(FCTDgrid.time,FCTDgrid.depth,FCTDgrid.temperature);
ax(1).CLim = [clim_temp(1) clim_temp(2)];
cb(1) = colorbar;
colormap(ax(1),cmocean('thermal'))
cb(1).Label.String = 'Temperature';

% Salinity
ax(2) = subtightplot(3,1,2);
pcolorjw(FCTDgrid.time,FCTDgrid.depth,real(FCTDgrid.salinity));
ax(2).CLim = [clim_sal(1) clim_sal(2)];
cb(2) = colorbar;
colormap(ax(2),cmocean('haline'))
cb(2).Label.String = 'Salinity';

% Chi
ax(3) = subtightplot(3,1,3);
pcolorjw(FCTDgrid.time,FCTDgrid.depth,log10(FCTDgrid.chi));
ax(3).CLim = [clim_chi(1) clim_chi(2)];
cb(3) = colorbar;
cb(3).Label.String = 'Chi';

% Add density contours and datetick
for iAx=1:3
   axes(ax(iAx))
   hold(ax(iAx),'on')
   [c,ch] = contour(FCTDgrid.time,FCTDgrid.depth,real(FCTDgrid.density),'k','levellist',levels_dens);
   contour(FCTDgrid.time,FCTDgrid.depth,FCTDgrid.temperature,'m','levellist',13);
   contour(FCTDgrid.time,FCTDgrid.depth,FCTDgrid.temperature,'c','levellist',15);
   clabel(c,ch);
   
   datetick(ax(iAx),'x','HH:MM','keeplimits')
end

% Depth axes
[ax(:).YLim] = deal([zlim(1) zlim(2)]);
[ax(:).YDir] = deal('reverse');

%MHA hack
%[ax(1:2).YLim] = deal([zlim(1) 500]);



