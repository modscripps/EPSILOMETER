function SD=mod_epsi_sd_buildtime(Meta_Data,a)

% a is the product of San's code for exemple see read_10_files_SODA_WW_d2.m)
% Meta_data is the the meta_data of the deployement
% 
% I am building a time axis in matlab time stamps
% we are using the EPSI timstamps (2st item in the EPSI block header)
% These time stamps resolution is sec. So I need to reconstruct the time
% axis. 
% Since the resolution is ~ 1 sec i often have twice the same stamps but
% not always
% 
% Issues arrive when we drop data block.
% I take the timestamp and go backward in time to build an array of 160
% timestamps
% 
% 

SD.madre=a.madre;
% epsi time

% get name channels. 
name_channels=Meta_Data.PROCESS.timeseries;

%get number of channel
nbchannels=Meta_Data.PROCESS.nb_channels;
% get start date
starttime=Meta_Data.starttime;
% get convert the MADRE timestamps in matlab stamps
timeheader=a.madre.TimeStamp/86400+datenum(1970,1,1);
timeheader=timeheader-timeheader(1)+starttime;

if isfield(SD.madre,'muTimeStamp')

    timeheader=[timeheader(1)-.5/86400; timeheader+SD.madre.muTimeStamp*30*30*1e-6/86400];
    SD.epsi.epsitime=nan(1,160*(numel(timeheader)-1));
    for t=1:numel(timeheader)-1
        SD.epsi.epsitime((t-1)*160+1:t*160)=linspace(timeheader(t)+1/320/86400,timeheader(t+1),160);
    end


    
else
    % time stamps on MADRE starts on january first 2017.
    %timeheader=starttime+(timeheader-datenum('01-01-2017 00:00:00'));
    %look for time difference between 0.5 second blocks
    dtimeheader=diff(timeheader);
    
    SD.epsi.epsitime=nan(1,160*numel(timeheader));
    
    last_t=0;
    count=0;
    flag_timebug=0;
    for t=1:numel(dtimeheader)
        dT=dtimeheader(t);
        if dT==0 || isnan(dT)
            count=count+1;
        end
        if dT>0
            if t==1
                SD.epsi.epsitime(1:160)=timeheader(1)-fliplr(linspace(1/325/86400,.5/86400,160));
            else
                if flag_timebug==0 % normal case
                    SD.epsi.epsitime((last_t+1)*160+1:(t+1)*160)=linspace(timeheader(last_t+1)+1/325/86400,timeheader(t+1),160+count*160);
                else   % if the timestamp bug and are decreasing
                    SD.epsi.epsitime((last_t+1)*160+1:(t+1)*160)=linspace(timeheader(t+1)-.5/86400,timeheader(t+1),160+count*160);
                end
            end
            last_t=t;
            count=0;
            flag_timebug=0;
        end
        if dT<0
            SD.epsi.epsitime((last_t+1)*160:(t+1)*160)=nan;
            last_t=t;
            count=0;
            flag_timebug=1;
            disp(t)
        end 
    end
    
    if dtimeheader(1)==0
        SD.epsi.epsitime(1:160)=SD.epsi.epsitime(161) - fliplr(linspace(1/325/86400,.5/86400,160));
    end
    if dtimeheader(end)==0
        SD.epsi.epsitime(end-(160+count*160)+1:end)=timeheader(end)- fliplr(linspace(1/325/86400,(.5+count*.5)/86400,160+count*160));
    end
    
    
end
SD.epsi.epsitime = SD.epsi.epsitime(:);

for n=1:nbchannels
    eval(sprintf('SD.epsi.%s=a.epsi.%s;',name_channels{n},name_channels{n}));
end

ind_OK=find(SD.epsi.epsitime>=SD.epsi.epsitime(1) & SD.epsi.epsitime<max(SD.epsi.epsitime));
SD.epsi.epsitime=SD.epsi.epsitime(ind_OK);
SD.epsi.EPSInbsample=a.epsi.EPSInbsample(ind_OK);
%epsitime=epsitime-epsitime(1);
for n=1:nbchannels
    eval(sprintf('SD.epsi.%s=SD.epsi.%s(ind_OK);',name_channels{n},name_channels{n}));
end
SD.epsi.flagSDSTR=SD.epsi.epsitime*0;

% commented to make GRANITE epsifish5sep09dep02
% can not really figure out if it will/wont work with the other
% deplyoments.
%SD.epsi.epsitime=linspace(SD.epsi.epsitime(1),SD.epsi.epsitime(end),length(SD.epsi.epsitime));

if isfield(Meta_Data,'SBEcal')
    % aux1 time
    % find the aux samples that matches the epsinbsample
    [~,iepsi1,iaux1] = intersect(a.epsi.EPSInbsample(ind_OK),a.aux1.Aux1Stamp);

    SD.aux1.T=a.aux1.T(iaux1);
    SD.aux1.P=a.aux1.P(iaux1);
    SD.aux1.C=a.aux1.C(iaux1);
    SD.aux1.S=sw_salt(SD.aux1.C*10./sw_c3515,SD.aux1.T,SD.aux1.P);
    SD.aux1.sig=sw_pden(SD.aux1.S,SD.aux1.T,SD.aux1.P,0);
    SD.aux1.aux1time=SD.epsi.epsitime(iepsi1);
end






