function [PT] = epsiProcess_get_profiles_from_PressureTimeseries(PressureTimeseries,Meta_Data,varargin)
% [PT] = epsiProcess_get_profiles_from_PressureTimeseries(PressureTimeseries,varargin)
%
% INPUTS
%   PressureTimeseries.dnum
%                     .time_s
%                     .P
%
% OPTIONAL INPUTS
% 'downLim','threshold','up','down','both'

plotFig = 0;
PT = PressureTimeseries;
p = PT.P;
% use nPoints median filter to smooth out the pressure field
%nPoints = Meta_Data.PROCESS.profs_nPoints;
%p = medfilt1(PT.P,nPoints);
% try to smooth out the data a bit
% dp = conv2(diff(conv2(p,ones(nPoints,1)/nPoints,'same'),1,1)',ones(1,nPoints)/nPoints,'same');
%
% NC - instead of first-differencing every point, first-difference over 5
% seconds. This way you can be more sure that the a profile is starting or
% ending. If you do that, you don't need to smooth the velocity data over
% much time. 32 points = 2 seconds, which can be significant at the bottom
% of the profile. Let's say nPoints = 16 for now.
numSec = 2;
nPoints = 16;
n_5sec = round(numSec/seconds(days(mode(diff(PT.dnum)))));
p_smoothed = conv2(p,ones(nPoints,1)/nPoints,'same');
for ii = 1:size(p_smoothed,1)-n_5sec
    dp_mid(ii) = [p_smoothed(ii+n_5sec) - p_smoothed(ii)]./numSec; %delta_P/delta_t where delta_t = numSec seconds
end
dp = interp1(linspace(0,1,length(dp_mid)),dp_mid,linspace(0,1,length(p_smoothed)));

% Set defaults (these will change in the next step if you specified in call
% to this function)
%downLim = 0.1;
%downLim = 0.025; (= 0.4m/s / 16 Hz). You want to choose something a little
%smaller than the fall rate, because you're detecting when the profile
%starts and it will take a short time to get up to speed.
downLim = Meta_Data.PROCESS.profs_downLim;
downLim = 0.1;
downCast = true;
minLength = Meta_Data.PROCESS.profs_minLength;

persistent argsNameToCheck;
if isempty(argsNameToCheck)
    argsNameToCheck = {'downLim','threshold','up','down','both'};
end
%
% index = 1;
% n_items = nargin-2;
% while (n_items > 0)
%     argsMatch = strcmpi(varargin{index},argsNameToCheck);
%     i = find(argsMatch,1);
%     if isempty(i)
%         error('MATLAB:epsiProcess_get_profiles_from_PressureTimeseries:wrongOption','Incorrect option specified: %s',varargin{index});
%     end
%
%     switch i
%         case 1 % downLim
%             if n_items == 1
%                 error('MATLAB:epsiProcess_get_profiles_from_PressureTimeseries:missingArgs','Missing input arguments');
%             end
%             downLim = varargin{index+1};
%             if downLim == 0
%                 error('MATLAB:epsiProcess_get_profiles_from_PressureTimeseries:downLim0Err','The threshold cannot be zero!');
%             end
%             index = index +2;
%             n_items = n_items-2;
%         case 2 % downLim
%             if n_items == 1
%                 error('MATLAB:epsiProcess_get_profiles_from_PressureTimeseries:missingArgs','Missing input arguments');
%             end
%             downLim = varargin{index+1};
%             if downLim == 0
%                 error('MATLAB:epsiProcess_get_profiles_from_PressureTimeseries:downLim0Err','The threshold cannot be zero!');
%             end
%             index = index +2;
%             n_items = n_items-2;
%         case 3 % up
%             downCast = false;
%
%             index = index + 1;
%             n_items = n_items - 1;
%         case 4 % down
%             downCast = true;
%
%             index = index + 1;
%             n_items = n_items - 1;
%     end
% end

%defining the threshold for going up and down
if downCast && downLim < 0
    downLim = -downLim;
    % elseif ~downCast && downLim > 0% going up
    %     downLim = -downLim;
end

if downCast
    dn = find(dp>downLim);
else
    dn = find(dp<downLim);
end

% Plot the criteria for dp
if plotFig
    figure

    % Plot pressure timeseries
    ax(1) = subplot(2,1,1);
    plot(PressureTimeseries.dnum,PressureTimeseries.P,'k','displayname','P');
    hold on
    plot(PressureTimeseries.dnum,p_smoothed,'b','displayname','smoothed P');
    set(gca,'ydir','reverse')
    title('P')
    datetick('x','HH:MM','keeplimits')

    % Plot dp
    dnum = PressureTimeseries.dnum;
    ax(2) = subplot(2,1,2);
    plot(dnum,dp,'b')
    title(sprintf('dp/d\_nSec, p smoothed over %3.0f points',nPoints))
    xx = get(gca,'xlim');
    hold on
    dl = plot([xx(1),xx(2)],[downLim downLim],'m');
    % NC - 'yline' is newer than Matlab 2018a
    %dl = yline(downLim,'m');
    legend(dl,sprintf('downLim = %3.3f', downLim),'location','best')
    datetick('x','HH:MM','keeplimits')
end



% find all indices of going down
if isempty(dn)
    return;
end
dn = [0, dn];

% find jumps in indices to indicate a start of a profile
startdown = dn(find(diff(dn)>1)+1);


if plotFig
    ax(1).NextPlot = 'add';
    plot(ax(1),PressureTimeseries.dnum(startdown),PressureTimeseries.P(startdown),'go','displayname','start downcast')

    ax(2).NextPlot = 'add';
    plot(ax(2),dnum(startdown),dp(startdown),'go','displayname','start downcast')

    linkprop([ax(:)],'xlim')
end

if isempty(startdown)
    return;
end

dn = dn(2:end);
try
    PT.drop = 0*PT.dnum;
catch
    PT.drop = 0*PT.time_s;
end

if dn(1)<startdown(1)
    startdown=[dn(1) startdown];
end

if startdown(end)<dn(end)
    startdown = [startdown dn(end)];
end
startprof = startdown + n_5sec; %Adjust startdown

% PT.drop is an array (the same size as other PT structure arrays)
% containing profile numbers.
enddown = nan(size(startdown,1),size(startdown,2));
for i=1:(length(startdown)-1)
    in = intersect(startdown(i):startdown(i+1)-1,dn);
    % All the up profiles are 0. All the down profiles are 1.
    PT.drop(in) = 1;
    enddown(i) = in(end);
end
enddown(i+1)=dn(end);
% Adjust enddown
endprof = enddown + round(n_5sec/2);


if plotFig
    ax(1).NextPlot = 'add';
    plot(ax(1),PressureTimeseries.dnum(endprof),PressureTimeseries.P(endprof),'ro','displayname','end downcast')

    ax(2).NextPlot = 'add';
    plot(ax(2),dnum(endprof),dp(endprof),'ro','displayname','end downcast')

    linkprop([ax(:)],'xlim')
end


% NC - Make sure profiles are at least a certain length
PT.startprof = startprof;
PT.endprof = endprof;
profLength = PT.P(PT.endprof)-PT.P(PT.startprof);
longEnough = abs(profLength)>=minLength;
PT.startprof = PT.startprof(longEnough);
PT.endprof = PT.endprof(longEnough);


if plotFig
    ax(1).NextPlot = 'add';
    plot(ax(1),PressureTimeseries.dnum(startprof(~longEnough)),PressureTimeseries.P(startprof(~longEnough)),'xk','displayname','too short')
    plot(ax(1),PressureTimeseries.dnum(endprof(~longEnough)),PressureTimeseries.P(endprof(~longEnough)),'xk','displayname','too short')

    ax(2).NextPlot = 'add';
    plot(ax(2),dnum(startprof(~longEnough)),dp(startprof(~longEnough)),'xk','displayname','too short')
    plot(ax(2),dnum(endprof(~longEnough)),dp(endprof(~longEnough)),'xk','displayname','too short')

    linkprop([ax(:)],'xlim')
end



if plotFig
    figure
    ax(1) = subplot(2,1,1);
    plot(PressureTimeseries.dnum,p_smoothed,'b','displayname','smoothed P');
    ax(1).NextPlot = 'add';

    plot(ax(1),PressureTimeseries.dnum(startprof(longEnough)),PressureTimeseries.P(startprof(longEnough)),'xg','displayname','too short')
    plot(ax(1),PressureTimeseries.dnum(endprof(longEnough)),PressureTimeseries.P(endprof(longEnough)),'xr','displayname','too short')

    ax(2) = subplot(2,1,2);
    plot(dnum,dp,'b')
    ax(2).NextPlot = 'add';

    plot(ax(2),dnum(startprof(longEnough)),dp(startprof(longEnough)),'xg','displayname','too short')

    plot(ax(2),dnum(endprof(longEnough)),dp(endprof(longEnough)),'xr','displayname','too short')

    linkprop([ax(:)],'xlim')
end



