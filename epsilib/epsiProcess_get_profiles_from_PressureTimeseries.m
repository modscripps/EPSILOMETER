function [PT] = epsiProcess_get_profiles_from_PressureTimeseries(PressureTimeseries,varargin)
 
% INPUTS
%   PressureTimeseries.dnum
%                     .time_s
%                     .P
%

plotFig = 0;
PT = PressureTimeseries;

% use nPoints median filter to smooth out the pressure field
nPoints = 64; % To do: make this an input. This would change if the CTD sampling rate changes
p = medfilt1(PT.P,nPoints);
% try to smooth out the data a bit
dp = conv2(diff(conv2(p,ones(nPoints,1)/nPoints,'same'),1,1)',ones(1,nPoints)/nPoints,'same');

%downLim = 0.1;
%downLim = 0.025;
downLim = 0.025; %minimum rate of change
downCast = true;
minLength = 20; % NC - Criteria for profile length - To do: add as an input

persistent argsNameToCheck;
if isempty(argsNameToCheck);
    argsNameToCheck = {'downLim','threshold','up','down','both'};
end

index = 1;
n_items = nargin-1;
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


% find all indices of going down
if isempty(dn)
    return;
end
dn = [0, dn];

% find jumps in indices to indicate a start of a profile
startdown = dn(find(diff(dn)>1)+1);

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


% Figure to check choices
if plotFig
    figure
    
    ax(1) = subtightplot(2,1,1);
    plot(PT.P)
    hold on
    plot(startdown,PT.P(PT.startprof),'*')
    plot(enddown,PT.P(PT.endprof),'o')
    set(gca,'ydir','reverse')
    
    ax(2) = subtightplot(2,1,2);
    plot(dp)
    hold on
    plot(dn,dp(dn))
    legend('dp','dp(dn)')
    
    linkprop([ax(:)],'xlim')
    a = [];
end


