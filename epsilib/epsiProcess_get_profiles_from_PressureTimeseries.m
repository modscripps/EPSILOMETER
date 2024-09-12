function [PT] = epsiProcess_get_profiles_from_PressureTimeseries(PressureTimeseries,Meta_Data)
% [PT] = epsiProcess_get_profiles_from_PressureTimeseries(PressureTimeseries,Meta_Data)
%
% INPUTS
%   PressureTimeseries.dnum
%                     .time_s
%                     .P
% Default processing values:
%   Meta_Data.PROCESS.numSec_lowpass = 4;
%                    .speedLim_down_start_m_s = 0.3;
%                    .speedLim_down_end_m_s = 0.1;
%                    .speedLim_up_start_m_s = 0.1;
%                    .speedLim_up_end_m_s = 0.05;
%                    .minLength_m = 10;
%                    .plotFig = 0;
% (All speed limits should be positive values in m/s (technically db/s) )

PT = PressureTimeseries;
p = PT.P;
dt = seconds(days(mode(diff(PT.dnum))));
sampling_rate_Hz = 1/dt;

% Set defaults
numSec_lowpass = 4;
% All speed limits should be positive values in m/s (technically db/s)
speedLim_down_start_m_s = 0.3;
speedLim_down_end_m_s = 0.1;
speedLim_up_start_m_s = 0.1;
speedLim_up_end_m_s = 0.05;
minLength_m = 10;
plotFig = 0;

% If Meta_Data has any of the profile processing fields specified, switch
% default values to those.
if isfield(Meta_Data,'PROFILES')
    if isfield(Meta_Data.PROFILES,'numSec_lowpass')
        numSec_lowpass = Meta_Data.PROFILES.numSec_lowpass;
    end
    if isfield(Meta_Data.PROFILES,'speedLim_down_start_m_s')
        speedLim_down_start_m_s = Meta_Data.PROFILES.speedLim_down_start_m_s;
    end
    if isfield(Meta_Data.PROFILES,'speedLim_down_end_m_s')
        speedLim_down_end_m_s = Meta_Data.PROFILES.speedLim_down_end_m_s;
    end
        if isfield(Meta_Data.PROFILES,'speedLim_up_start_m_s')
        speedLim_up_start_m_s = Meta_Data.PROFILES.speedLim_up_start_m_s;
    end
    if isfield(Meta_Data.PROFILES,'speedLim_up_end_m_s')
        speedLim_up_end_m_s = Meta_Data.PROFILES.speedLim_up_end_m_s;
    end
    if isfield(Meta_Data.PROFILES,'minLength_m')
        minLength_m = Meta_Data.PROFILES.minLength_m;
    end
    if isfield(Meta_Data.PROFILES,'plotFig')
        plotFig = Meta_Data.PROFILES.plotFig;
    end
end

% % Convert number of seconds to number of samples
% n_smooth = round(numSec_smooth*sampling_rate_Hz);
% n_1stdiff = round(numSec_1stdiff*sampling_rate_Hz);

% Find freuencies for applying lowpass filter
lowpass_cutoff_Hz = 1/numSec_lowpass;
Fc=lowpass_cutoff_Hz./(sampling_rate_Hz/2);
[b,a]=cheby2(3,20,Fc,"low");

% Lowpass filter the pressure, and get dpdt
% p_lowpass = lowpass(fillgaps(p),lowpass_cutoff_Hz,sampling_rate_Hz);
p_lowpass = filtfilt(b,a,fillgaps(p));
dpdt_mid = diff(p_lowpass)./dt;
dpdt = interp1(linspace(0,1,length(dpdt_mid)),dpdt_mid,linspace(0,1,length(p)));


% % Smooth pressure and get dpdt
% p_smoothed = movmean(p,n_smooth);
% for ii = 1:size(p_smoothed,1)-n_1stdiff
%     dpdt_mid(ii) = [p_smoothed(ii+n_1stdiff) - p_smoothed(ii)]./numSec_1stdiff;
% end
% dpdt = interp1(linspace(0,1,length(dpdt_mid)),dpdt_mid,linspace(0,1,length(p_smoothed)));
%
% % Lowpass filter dpdt

% Initialize a 'drop' field in PressureTimeseries for counting profiles.
try
    PT.drop = 0*PT.dnum;
catch
    PT.drop = 0*PT.time_s;
end

% First get downcasts, then get upcasts
for downCast=[1,0]
    clear profiling_start profiling_end startcast endcast

    % % Determine the appropriate speed limits based on the direction
    % if downCast
    %     DPDT = dpdt;
    %     start_speed_lim = speedLim_down_start_m_s;
    %     end_speed_lim = speedLim_down_end_m_s;
    % else
    %     DPDT = -dpdt;
    %     start_speed_lim = speedLim_up_start_m_s;
    %     end_speed_lim = speedLim_up_end_m_s;
    % end


    % Find indices for the start and end of downcasts
    if downCast
        profiling_start_criteria = find(dpdt>0 & dpdt>speedLim_down_start_m_s);
        profiling_end_criteria = find(dpdt>0 & dpdt>speedLim_down_end_m_s);
    else %upcast
        % Find indices for the start and end of upcasts
        profiling_start_criteria = find(-dpdt>0 & -dpdt>speedLim_up_start_m_s);
        profiling_end_criteria = find(-dpdt>0 & -dpdt>speedLim_up_end_m_s);
    end


    % % Find where pressure differences are above the starting speed limit
    % profiling_start_idx = DPDT > start_speed_lim;
    % profiling_start = find(profiling_start_idx);
    % % Find where pressure differences are above the ending speed limit, but
    % % ignore the indices that meet the starting speed limit
    % profiling_end_idx = DPDT < end_speed_lim;
    % profiling_end = find(profiling_end_idx);

    % Find jumps in indices and add 1 to indicate a start of a profile
    keep = [1, find(diff([profiling_start_criteria,numel(profiling_start_criteria)+10]) > 1) + 1];
    % 2024 05 18 ASTRAL hacking to make it work, Need reevaluation
    try
    profiling_start = profiling_start_criteria(keep);
    catch
        profiling_start = [];
    end

    % Find jumps in indiced to indicate end of profile
    keep = [find(diff([profiling_end_criteria,numel(profiling_start_criteria)+10]) > 1),length(profiling_end_criteria)];
    try
        profiling_end = sort(profiling_end_criteria(keep));
    catch
        profiling_end = [];
    end

    % Initialize variables to store start and end points of profiles
    startcast = [];
    endcast = [];

    % Iterate through each start point to find the corresponding end point
    for i = 1:length(profiling_start)
        % Find the first end point that occurs after the start point
        end_idx = find(profiling_end > profiling_start(i), 1);

        % If an end point exists, store the start and end points
        if ~isempty(end_idx)
            startcast = [startcast; profiling_start(i)];
            endcast = [endcast; profiling_end(end_idx)];
        end
    end

    % It's possible that the same profiling_end is chosen for multiple
    % profiling_starts. unique() will choose the last profiling_start.
    % I'm not sure if it's better to keep the first or last, but let's go
    % with the last because it's easy to get
    [~,iU] = unique(endcast);
    startcast = startcast(iU);
    endcast = endcast(iU);

    % Make sure profiles are at least a certain length
    profLength = abs(p(endcast) - p(startcast));
    longEnough = profLength >= minLength_m;

    % Store the start and end indices of profiles based on direction
    if downCast
        PT.startdown = startcast(longEnough);
        PT.enddown = endcast(longEnough);

    else
        PT.startup = startcast(longEnough);
        PT.endup = endcast(longEnough);
    end

    if plotFig
        ax = make_cast_plots_1(downCast,PT,dpdt,profiling_start_criteria,profiling_start,profiling_end_criteria,profiling_end);
        % Link axes so you can zoom in
        if downCast
            d1 = linkprop([ax(:)],'xlim');
            d2 = linkprop([ax([1,4])],'ylim');
            d3 = linkprop([ax([2,5])],'ylim');
            d4 = linkprop([ax([3,6])],'ylim');
        else
            u1 = linkprop([ax(:)],'xlim');
            u2 = linkprop([ax([1,4])],'ylim');
            u3 = linkprop([ax([2,5])],'ylim');
            u4 = linkprop([ax([3,6])],'ylim');
        end
    end

end %End loop through downcasts then upcasts

% Make sure the profiles are always alternating down, up, down, up. If
% there are two downs in a row or two ups in a row, merge those profiles.
%
% TODO: This assumes the profiler goes down first and then up. If it goes
% up first, swap the criteria.
%
% Merge all upcast and downcast indices - you'll use these to look for up-
% or downcasts in between potential gaps
%
% Do this a few times in case there are more than one successive up or down
% casts
repeat_i = 1;
while repeat_i < 10
    all_up = [];
    all_down = [];
    for ii=1:length(PT.startup)
        all_up = [all_up, PT.startup(ii):PT.endup(ii)];
    end
    for ii=1:length(PT.startdown)
        all_down = [all_down, PT.startdown(ii):PT.enddown(ii)];
    end

    % Look for upcasts in a row. Merge if they are
    ii = 2;
    while ii<numel(PT.startup)
        if ~any(intersect(PT.endup(ii-1):PT.startup(ii),all_down))
            % Merge upcasts ii and ii-1
            PT.startup(ii-1) = PT.startup(ii-1);
            PT.endup(ii-1) = PT.endup(ii);
            PT.startup(ii) = [];
            PT.endup(ii) = [];
            ii = ii+1;
        else
            ii = ii+1;
        end
    end

    % Look for downcasts in a row. Merge if they are
    ii = 2;
    while ii<numel(PT.startdown)
        if ~any(intersect(PT.enddown(ii-1):PT.startdown(ii),all_up))
            % Merge downcasts ii and ii-1
            PT.startdown(ii-1) = PT.startdown(ii-1);
            PT.enddown(ii-1) = PT.enddown(ii);
            PT.startdown(ii) = [];
            PT.enddown(ii) = [];
            ii = ii+1;
        else
            ii = ii+1;
        end
    end
    repeat_i = repeat_i + 1;
end %end repeat_i


%% Populate PT.drop
for i=1:length(PT.startdown)
    in = PT.startdown(i):PT.enddown(i);
    PT.drop(in) = i;
end
for i=1:length(PT.startup)
    in = PT.startup(i):PT.endup(i);
    PT.drop(in) = -i;
end

%% If you're using downcasts startprof/endprof = startdown/enddown. 
% If you're using upcasts, startprof/endprof = startup/endup
% If you're using both, startprof/endprof = startdown/enddown and startup/endup
switch Meta_Data.paths.profiles
    case 'down'
        PT.startprof = PT.startdown;
        PT.endprof = PT.enddown;
    case 'up'
        PT.startprof = PT.startup;
        PT.endprof = PT.endup;
    case 'both'
        PT.startprof = sort([PT.startdown(:);PT.startup(:)]);
        PT.endprof = sort([PT.enddown(:);PT.endup(:)]);
end


%% Plot the criteria for dp and selected profiles
if plotFig

    ax = make_cast_plots_2(PT,p_lowpass,dpdt,speedLim_down_start_m_s,speedLim_down_end_m_s,speedLim_up_start_m_s,speedLim_up_end_m_s);
    lp = linkprop([ax(:)],'xlim');

end

end %/end epsiProcess_get_profiles_from_PressureTimeseries


%% SUBFUNCTIONS FOR MAKING PLOTS
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------


function [ax] = make_cast_plots_1(downCast,PT,dpdt,profiling_start_criteria,profiling_start,profiling_end_criteria,profiling_end)
clear paired
get_paired

figure('units','inches','position',[0 0 23 13])

idx = 1:length(PT.dnum);

% ax(1) = subtightplot(4,1,1);
% p1(1) = plot(PT.dnum(profiling_start_criteria),idx(profiling_start_criteria),'o','color',paired.lightPurple);
% hold on
% p1(2) = plot(PT.dnum(profiling_start),idx(profiling_start),'x','color',paired.purple);
% grid on
% ylabel('index #')
%

% START CASTS
ax(1) = subtightplot(3,2,1);
p2(1) = plot(PT.dnum(profiling_start_criteria),diff([profiling_start_criteria,0]),'o','color',paired.lightPurple);
%ylim([-32 60*16])
%ytick(0:16*5:60*16)
grid on
ylabel('diff(index #)')
if downCast
    title('Start of downcasts')
else
    title('Start of upcasts')
end

ax(2) = subtightplot(3,2,3);
plot(PT.dnum, dpdt,'.','color',[0.5 0.5 0.5])
hold on
plot(PT.dnum(profiling_start_criteria), dpdt(profiling_start_criteria),'o','color',paired.lightPurple);
plot(PT.dnum(profiling_start), dpdt(profiling_start),'x','color',paired.purple)
if downCast
    plot(PT.dnum(PT.startdown), dpdt(PT.startdown),'o','color','m')
else
    plot(PT.dnum(PT.startup), dpdt(PT.startup),'o','color','m')
end
ylabel('dpdt');

ax(3) = subtightplot(3,2,5);
plot(PT.dnum, PT.P,'.','color',[0.5 0.5 0.5])
hold on
plot(PT.dnum(profiling_start_criteria), PT.P(profiling_start_criteria),'o','color',paired.lightPurple)
plot(PT.dnum(profiling_start), PT.P(profiling_start),'x','color',paired.purple)
if downCast
    plot(PT.dnum(PT.startdown), PT.P(PT.startdown),'o','color','m')
else
    plot(PT.dnum(PT.startup), PT.P(PT.startup),'o','color','m')
end
ylabel('pressure')
ax(2).YDir = 'reverse';


% ENDCASTS
ax(4) = subtightplot(3,2,2);
p2(1) = plot(PT.dnum(profiling_end_criteria),diff([profiling_end_criteria,0]),'o','color',paired.lightGreen);
%ylim([-32 60*16])
%ytick(0:16*5:60*16)
grid on
if downCast
    title('End of downcasts')
else
    title('End of upcasts')
end

ax(5) = subtightplot(3,2,4);
plot(PT.dnum, dpdt,'.','color',[0.5 0.5 0.5])
hold on
plot(PT.dnum(profiling_end_criteria), dpdt(profiling_end_criteria),'o','color',paired.lightGreen);
plot(PT.dnum(profiling_end), dpdt(profiling_end),'x','color',paired.green)
if downCast
    plot(PT.dnum(PT.enddown), dpdt(PT.enddown),'o','color','c')
else
    plot(PT.dnum(PT.endup), dpdt(PT.endup),'o','color','c')
end
grid on

ax(6) = subtightplot(3,2,6);
plot(PT.dnum, PT.P,'.','color',[0.5 0.5 0.5])
hold on
plot(PT.dnum(profiling_end_criteria), PT.P(profiling_end_criteria),'o','color',paired.lightGreen)
plot(PT.dnum(profiling_end), PT.P(profiling_end),'x','color',paired.green)
if downCast
    plot(PT.dnum(PT.enddown), PT.P(PT.enddown),'o','color','c')
else
    plot(PT.dnum(PT.endup), PT.P(PT.endup),'o','color','c')
end

[ax([3,6]).YDir] = deal('reverse');
[ax(2:6).XLim] = deal(ax(1).XLim);
[ax(:).XTick] = deal('');
for iAx=5:6
    datetick(ax(iAx),'x','HH:MM','keeplimits');
end

lp1 = linkprop([ax(:)],'xlim');
lp3 = linkprop([ax([1,4])],'ylim');
lp3 = linkprop([ax([2,5])],'ylim');
lp4 = linkprop([ax([3,6])],'ylim');
end %/end make_cast_plots_1

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

function [ax] = make_cast_plots_2(PT,p_lowpass,dpdt,speedLim_down_start_m_s,speedLim_down_end_m_s,speedLim_up_start_m_s,speedLim_up_end_m_s)

% Get colors for pairs of profiles - stack 10 'paired' colors 500 times
paired10 = colorbrewer('Paired');
paired10 = paired10(1:10,:);
paired = repmat(paired10,500,1);
downColors = paired(2:2:end,:);
upColors = paired(1:2:end,:);

%
figure('units','inches','position',[0 0 23 13])
gap = [0.05 0.05];
marg_v = [0.05 0.05];
marg_h = [0.05 0.05];

ax(1) = subtightplot(2,2,1,gap,marg_v,marg_h);
plot(PT.dnum,PT.P,'k.');
hold on
plot(PT.dnum,p_lowpass,'color',[0.5 0.5 0.5])
set(gca,'ydir','reverse');
hold on
for iP=1:max(abs(PT.drop))
    down = PT.drop==iP;
    up = PT.drop==-iP;
    s1 = scatter(PT.dnum(down),PT.P(down),'o','markeredgecolor',downColors(iP,:));
    s2 = scatter(PT.dnum(up),PT.P(up),'o','markeredgecolor',upColors(iP,:));
end
ylabel('Depth (m)')
try
    legend([s1,s2],{'downcasts','upcasts'})
catch
end

ax(2) = subtightplot(2,2,3,gap,marg_v,marg_h);
plot(PT.dnum,dpdt,'k.');
hold on
for iP=1:max(PT.drop)
    down = PT.drop==iP;
    up = PT.drop==-iP;
    scatter(PT.dnum(down),dpdt(down),'o','markeredgecolor',downColors(iP,:));
    scatter(PT.dnum(up),dpdt(up),'o','markeredgecolor',upColors(iP,:));
end
ln.down_start = yline(speedLim_down_start_m_s,'m--','linewidth',1,'displayname','down start');
ln.up_start = yline(-speedLim_up_start_m_s,'c--','linewidth',1,'displayname','up start');
ln.down_end = yline(speedLim_down_end_m_s,'m-','linewidth',1,'displayname','down end');
ln.up_end = yline(-speedLim_up_end_m_s,'c-','linewidth',1,'displayname','up end');
ylabel('')
legend([ln.down_start,ln.down_end,ln.up_start,ln.up_end])
ylabel('Speed (m/s)')

ax(3) = subtightplot(2,2,2,gap,marg_v,marg_h);
plot(PT.dnum,PT.P,'k','displayname','P');
hold on
plot(PT.dnum,p_lowpass,'color',[0.5 0.5 0.5],'displayname','smoothed P');
ax(3).NextPlot = 'add';
plot(ax(3),PT.dnum(PT.startdown),PT.P(PT.startdown),'om','displayname','start down')
plot(ax(3),PT.dnum(PT.enddown),PT.P(PT.enddown),'xm','displayname','end down')
plot(ax(3),PT.dnum(PT.startup),PT.P(PT.startup),'oc','displayname','start up')
plot(ax(3),PT.dnum(PT.endup),PT.P(PT.endup),'xc','displayname','end up')

ax(4) = subtightplot(2,2,4,gap,marg_v,marg_h);
plot(PT.dnum,dpdt,'k')
ax(4).NextPlot = 'add';
plot(ax(4),PT.dnum(PT.startdown),dpdt(PT.startdown),'om')
plot(ax(4),PT.dnum(PT.enddown),dpdt(PT.enddown),'xm')
plot(ax(4),PT.dnum(PT.startup),dpdt(PT.startup),'oc')
plot(ax(4),PT.dnum(PT.endup),dpdt(PT.endup),'xc')
ln.down = yline(speedLim_down_start_m_s,'m--','linewidth',1,'displayname','down start');
ln.up = yline(-speedLim_up_start_m_s,'c--','linewidth',1,'displayname','up start');
ln.down = yline(speedLim_down_end_m_s,'m-','linewidth',1,'displayname','down end');
ln.up = yline(-speedLim_up_end_m_s,'c-','linewidth',1,'displayname','up end');

for iAx=1:4
    datetick(ax(iAx),'x','keeplimits');
end
[ax([1,3]).YDir] = deal('reverse');
[ax([2,4]).YLim] = deal([-1 1]);

linkprop([ax(:)],'xlim')

try
    figname = fullfile(Meta_Data.paths.figures,'pick_profiles');
    savefig(figname)
catch
    disp('Failed to save pick_profiles.fig - epsiProcess_get_profiles_from_PressureTimeseries.m');
end

end %/end make_cast_plots_2
