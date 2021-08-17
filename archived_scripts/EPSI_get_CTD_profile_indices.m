function [CTD] = EPSI_get_CTD_profile_indices(MATdir)
% New version of EPSI_create_profiles  that finds profiles based on dPdz
% anywhere the fish turns around, instead of requiring it to come to the
% surface. Find profiles using the ctd structure and interpolate the epsi
% data onto those profiles. 
% 
% TODO:
%   get ups and downs too and choose profiles based on whether you want
%   ups, downs, or both
%
% INPUTS:
%   ctd structure
%
% OPTIONAL INPUT FLAGS
%   {'downLim','threshold','upcast','downcast'};
% 
% Edited from:
% FCTD = FastCTD_FindCasts(FCTD);
%   finds when the FastCTD is going up and when it is going down
%   if 'UPCAST' is specified then the profile will search for upcast but by
%   defalt it would search for downcasts
%
%   In this data set we assume the data is in descending order of time
%
% Written by Jody Klymak
% Updated 2011 07 14 by San Nguyen
% Updated 2012 09 29 by San Nguyen

plotFig = 0;

% Load PressureTimeseries structure from MATdir
load(fullfile(MATdir,'Epsi_PressureTimeseries'))

CTD = PressureTimeseries;

if ~isfield(CTD,'P')
    disp('There is no PressureTimeseries.P')
    return;
end

% use nPoints median filter to smooth out the pressure field
nPoints = 64; 
p = medfilt1(CTD.P,nPoints);
% try to smooth out the data a bit
dp = conv2(diff(conv2(p,ones(nPoints,1)/nPoints,'same'),1,1)',ones(1,nPoints)/nPoints,'same');

%downLim = 0.1;
%downLim = 0.025;
downLim = 0.025;
downCast = true;
minLength = 20; % NC - Criteria for profile length - To do: add as an input

persistent argsNameToCheck;
if isempty(argsNameToCheck);
    argsNameToCheck = {'downLim','threshold','upcast','downcast'};
end

index = 1;
n_items = nargin-1;
while (n_items > 0)
    argsMatch = strcmpi(varargin{index},argsNameToCheck);
    i = find(argsMatch,1);
    if isempty(i)
        error('MATLAB:FastCTD_FindCasts:wrongOption','Incorrect option specified: %s',varargin{index});
    end
    
    switch i
        case 1 % downLim
            if n_items == 1
                error('MATLAB:FastCTD_FindCasts:missingArgs','Missing input arguments');
            end
            downLim = varargin{index+1};
            if downLim == 0 
                error('MATLAB:FastCTD_FindCasts:downLim0Err','The threshold cannot be zero!');
            end
            index = index +2;
            n_items = n_items-2;
        case 2 % downLim
            if n_items == 1
                error('MATLAB:FastCTD_FindCasts:missingArgs','Missing input arguments');
            end
            downLim = varargin{index+1};
            if downLim == 0 
                error('MATLAB:FastCTD_FindCasts:downLim0Err','The threshold cannot be zero!');
            end
            index = index +2;
            n_items = n_items-2;
        case 3 % upcast
            downCast = false;
            
            index = index + 1;
            n_items = n_items - 1;
        case 4 % downcast
            downCast = true;
            
            index = index + 1;
            n_items = n_items - 1;
    end
end

%defining the threshold for going up and down
if downCast && downLim < 0
    downLim = -downLim;
elseif (~downCast) && downLim > 0% going up
    downLim = -downLim;
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
    CTD.drop = 0*CTD.ctddnum;
catch
    CTD.drop = 0*CTD.ctdtime;
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


% CTD.drop is an array (the same size as other CTD strucutre arrays)
% containing profile numbers.
enddown = nan(size(startdown,1),size(startdown,2));
for i=1:(length(startdown)-1)
    in = intersect(startdown(i):startdown(i+1)-1,dn);  
    % All the up profiles are 0. All the down profiles are 1.
    CTD.drop(in) = 1;
    enddown(i) = in(end);
end
enddown(i+1)=dn(end);


% NC - Make sure profiles are at least a certain length
CTD.startprof = startdown;
CTD.endprof = enddown;
profLength = CTD.P(CTD.endprof)-CTD.P(CTD.startprof);
longEnough = profLength>=minLength;
CTD.startprof = CTD.startprof(longEnough);
CTD.endprof = CTD.endprof(longEnough);

% Figure to check choices
if plotFig
    figure
    
    ax(1) = subtightplot(2,1,1);
    plot(CTD.P)
    hold on
    plot(startdown,CTD.P(CTD.startprof),'*')
    plot(enddown,CTD.P(CTD.endprof),'o')
    set(gca,'ydir','reverse')
    
    ax(2) = subtightplot(2,1,2);
    plot(dp)
    hold on
    plot(dn,dp(dn))
    legend('dp','dp(dn)')
    
    linkprop([ax(:)],'xlim')
    a = [];
end



PressureTimeseries = CTD;
save(fullfile(MATdir,'Epsi_PressureTimeseries'),'PressureTimeseries')

end

