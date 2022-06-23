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

% use nPoints median filter to smooth out the pressure field
nPoints = Meta_Data.PROCESS.profs_nPoints;
p = medfilt1(PT.P,nPoints);
% try to smooth out the data a bit
dp = conv2(diff(conv2(p,ones(nPoints,1)/nPoints,'same'),1,1)',ones(1,nPoints)/nPoints,'same');

% Set defaults (these will change in the next step if you specified in call
% to this function)
%downLim = 0.1;
%downLim = 0.025; (= 0.4m/s / 16 Hz). You want to choose something a little
%smaller than the fall rate, because you're detecting when the profile
%starts and it will take a short time to get up to speed.
downLim = Meta_Data.PROCESS.profs_downLim;
downCast = true;
minLength = Meta_Data.PROCESS.profs_minLength;

persistent argsNameToCheck;
if isempty(argsNameToCheck)
    argsNameToCheck = {'downLim','threshold','up','down','both'};
end

index = 1;
n_items = nargin-2;
while (n_items > 0)
    argsMatch = strcmpi(varargin{index},argsNameToCheck);
    i = find(argsMatch,1);
    if isempty(i)
        error('MATLAB:epsiProcess_get_profiles_from_PressureTimeseries:wrongOption','Incorrect option specified: %s',varargin{index});
    end
    
    switch i
        case 1 % downLim
            if n_items == 1
                error('MATLAB:epsiProcess_get_profiles_from_PressureTimeseries:missingArgs','Missing input arguments');
            end
            downLim = varargin{index+1};
            if downLim == 0 
                error('MATLAB:epsiProcess_get_profiles_from_PressureTimeseries:downLim0Err','The threshold cannot be zero!');
            end
            index = index +2;
            n_items = n_items-2;
        case 2 % downLim
            if n_items == 1
                error('MATLAB:epsiProcess_get_profiles_from_PressureTimeseries:missingArgs','Missing input arguments');
            end
            downLim = varargin{index+1};
            if downLim == 0 
                error('MATLAB:epsiProcess_get_profiles_from_PressureTimeseries:downLim0Err','The threshold cannot be zero!');
            end
            index = index +2;
            n_items = n_items-2;
        case 3 % up
            downCast = false;
            
            index = index + 1;
            n_items = n_items - 1;
        case 4 % down
            downCast = true;
            
            index = index + 1;
            n_items = n_items - 1;
    end
end

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
    % Plot pressure timeseries
    ax(1) = subplot(2,1,1);
    plot(PressureTimeseries.dnum,PressureTimeseries.P,'k');
    set(gca,'ydir','reverse')
    title('P')
    datetick('x','HH:MM','keeplimits')

    % Plot dp
    mid_dnum = nanmean([PressureTimeseries.dnum(1:end-1),...
        PressureTimeseries.dnum(2:end)],2);
    ax(2) = subplot(2,1,2);
    plot(mid_dnum,dp,'b')
    title(sprintf('dp, p smoothed over %3.0f points',nPoints))
    dl = yline(downLim,'m');
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
plot(ax(1),PressureTimeseries.dnum(startdown),PressureTimeseries.P(startdown),'go')

ax(2).NextPlot = 'add';
plot(ax(2),mid_dnum(startdown),dp(startdown),'go')
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

% % NC - Make sure startdowns are a certain distance apart
% distBetweenStartdowns = 1*60*16; %(1 minutes at 16Hz sampling)
% startdowndiff = [0,diff(startdown)];
% idxKeep = startdowndiff>distBetweenStartdowns;
% startdown = startdown(idxKeep);


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


% NC - Make sure profiles are at least a certain length
PT.startprof = startdown;
PT.endprof = enddown;
profLength = PT.P(PT.endprof)-PT.P(PT.startprof);
longEnough = abs(profLength)>=minLength;
PT.startprof = PT.startprof(longEnough);
PT.endprof = PT.endprof(longEnough);


