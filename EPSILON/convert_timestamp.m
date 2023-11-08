function [time_s,dnum] = convert_timestamp(timestamp)
seconds1970 = timestamp./1000;
days1970 = seconds1970/(24*60*60);
offset1970_days = datenum(1970,1,1) - datenum(0000,1,0);
offset1970_seconds = offset1970_days.*(24*60*60);
time_s = seconds1970 + offset1970_seconds;
dnum = days1970 + offset1970_days;
end