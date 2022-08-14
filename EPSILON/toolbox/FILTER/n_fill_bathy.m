function [f] = n_fill_bathy(x,y)       

% Make x and y row vectors
x = x(:).';
y = y(:).';

ax = gca;
ax.NextPlot = 'add';

xLim = ax.XLim;
yLim = ax.YLim;
% Technique 1: Expand out from the axes limits
xLow = xLim(1)-10;
yLow = yLim(2)+10;
% Technique 2: Expand out from the data limits
xLow = nanmin(x);
yLow = nanmax(y)+200;

% Find bathy
nans = isnan(x) | isnan(y);
x(nans) = [];
y(nans) = [];

fill([xLow,x,x(end),xLow,xLow],[y(1),y,yLow,yLow,y(1)],[0.5 0.5 0.5]);