function time_round=round_time_AAA(varargin)
% function round_time_AAA rounds a deisred time to a multiple of the selected unit
%
% Inputs:
%   time - matlab datenum time, if empty uses 'now'
%   unit - 'year','month','day','hour','minute','second'
%   value - time step in your selected unit
%   direction - 'round','floor','ceil'
%
% Outputs:
%   time_round - rounded time
%
% Example:
%   time_round = round_time_AAA(time,'minute',15,'floor'); rounds down to
%   the nearest 15 minutes starting on the hour
%
%   time_round = round_time_AAA(time,'second',10,'round'); rounds to the
%   nearest 10 seconds starting on the minute
%
%   time_round = round_time_AAA(time,'year',50,'round'); rounds to the
%   nearest 50 years starting at 0 C.E.
% 
%   time_round = round_time_AAA(time,'day',50,'round'); rounds to the
%   nearest 50 yeardays starting from january 1st
%
%   time_round = round_time_AAA; returns the current time to the nearest
%   second
%
% Notes:
%   This only works for a single value time input, to apply this to multiple times run in a for loop.
%
% Alex Andriatis
% 2023-11-10

P=inputParser;

time=now;
addOptional(P,'time',time,@isnumeric);

checkString=@(s) any(strcmp(s,{'year','month','day','hour','minute','second'}));
addOptional(P,'unit','second',checkString);

addOptional(P,'value',1,@isnumeric);

checkString=@(s) any(strcmp(s,{'round','floor','ceil'}));
addOptional(P,'direction','round',checkString);

parse(P,varargin{:});
time = P.Results.time;
unit = P.Results.unit;
value = P.Results.value;
direction = P.Results.direction;

switch unit
    case 'year'
        t1=year(time);
        if t1>0
            tgrid=datenum(0:value:t1+value,1,1);
        else
            tgrid=datenum(0:-value:t1-value,1,1);
            tgrid=fliplr(tgrid);
        end
    case 'month'
        t1=month(time);
        tgrid=datenum(year(time),1:value:t1+value,1);
    case 'day'
        t1=floor(time-datenum(year(time),1,1))+1;
        tgrid=datenum(year(time),1,1:value:t1+value);
    case 'hour'
        t1=hour(time);
        tgrid=datenum(year(time),month(time),day(time),0:value:t1+value,0,0);
    case 'minute'
        t1=minute(time);
        tgrid=datenum(year(time),month(time),day(time),hour(time),0:value:t1+value,0);
    case 'second'
        t1=second(time);
        tgrid=datenum(year(time),month(time),day(time),hour(time),minute(time),0:value:t1+value);
end

switch direction
    case 'round'
        [~,It]=min(abs(tgrid-time));
    case 'floor'
        It=find((tgrid-time)<=0,1,'last');    
    case 'ceil'
        It=find((tgrid-time)>=0,1,'first');
end
time_round=tgrid(It);