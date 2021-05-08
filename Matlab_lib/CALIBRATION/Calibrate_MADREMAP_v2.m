function Calibrate_MADREMAP_v2(Meta_Data,tscan)

%  script to calibrate the MADRE-MAP set. Good practice try to keep MADRE
%  and MAP together !!!!!
%  Should be use with data for 10 min tests in the Faraday cage with dummy
%  probes (TODO: define dummy probes)
% 

%  input: Meta_Data, file
%  created with Meta_Data=create_Meta_Data(file). Meta_Data contain the
%  path to calibration file and EPSI configuration needed to process the
%  epsi data
%
%  file is the file created after 
%  Created by Arnaud Le Boyer on 7/28/18.
%  Copyright © 2018 Arnaud Le Boyer. All rights reserved.




addpath ../toolbox/

%file='/Users/aleboyer/ARNAUD/SCRIPPS/EPSILOMETER/CALIBRATION/ELECTRONICS/MADRE_MAP/MM1/epsi/MM1.mat'
%file='/Volumes/DataDrive/GRANITE/For_Arnaud_Sept2018/benchtest_for_alb/aug30_11/epsi/epsi_aug30_11.mat'
%load('/Users/aleboyer/ARNAUD/SCRIPPS/EPSILOMETER/CALIBRATION/ELECTRONICS/MADRE_MAP/MM1/raw/Meta_DEV_MM1.mat')
EPSI=load(fullfile(Meta_Data.Epsipath,['epsi_' Meta_Data.deployment '.mat']));

xlabel_str=datestr(now,'mm-dd-yyyy');
timeaxis=86400*(EPSI.epsitime-EPSI.epsitime(1));
EPSI.ramp_count=EPSI.a1.*0;
dramp=diff(EPSI.ramp_count);

L=length(timeaxis);

close all
%figure('units','inch','position',[0,0,35,15]);
figure
cmap=colormap(parula(8));
figure(1)
ax(1)=subplot('Position',[.1 .9 .8 .05]);
ax(2)=subplot('Position',[.1 .84 .8 .05]);
ax(3)=subplot('Position',[.1 .78 .8 .05]);
ax(4)=subplot('Position',[.1 .72 .8 .05]);
ax(5)=subplot('Position',[.1 .66 .8 .05]);
plot(ax(1),timeaxis(floor(L/2)-325*50:floor(L/2)+325*50),EPSI.a1(floor(L/2)-325*50:floor(L/2)+325*50),'Color',cmap(1,:))
hold(ax(1),'on')
plot(ax(1),timeaxis(floor(L/2)-325*50:floor(L/2)+325*50),EPSI.a2(floor(L/2)-325*50:floor(L/2)+325*50),'Color',cmap(2,:))
hold(ax(1),'off')

plot(ax(2),timeaxis(floor(L/2)-325*50:floor(L/2)+325*50),EPSI.a3(floor(L/2)-325*50:floor(L/2)+325*50),'Color',cmap(3,:))

plot(ax(3),timeaxis(floor(L/2)-325*50:floor(L/2)+325*50),EPSI.t1(floor(L/2)-325*50:floor(L/2)+325*50),'Color',cmap(4,:))
hold(ax(3),'on')
plot(ax(3),timeaxis(floor(L/2)-325*50:floor(L/2)+325*50),EPSI.t2(floor(L/2)-325*50:floor(L/2)+325*50),'Color',cmap(5,:))
hold(ax(3),'off')

plot(ax(4),timeaxis(floor(L/2)-325*50:floor(L/2)+325*50),EPSI.s1(floor(L/2)-325*50:floor(L/2)+325*50),'Color',cmap(6,:))
hold(ax(4),'on')
plot(ax(4),timeaxis(floor(L/2)-325*50:floor(L/2)+325*50),EPSI.s2(floor(L/2)-325*50:floor(L/2)+325*50),'Color',cmap(7,:))
hold(ax(4),'off')

plot(ax(5),timeaxis(floor(L/2)-325*50:floor(L/2)+325*50),dramp(floor(L/2)-325*50:floor(L/2)+325*50),'Color',cmap(8,:))



 ax(1).YLim=[-.05 .52];
 ax(2).YLim=[-1.2 -.6];
 ax(3).YLim=[1.2 1.3];
 ax(4).YLim=[1.2 1.3];
 ax(5).YLim=[0 2];
ylabel(ax(1),'g','FontSize',20)
ylabel(ax(2),'g','FontSize',20)
ylabel(ax(3),'V','FontSize',20)
ylabel(ax(4),'V','FontSize',20)
ylabel(ax(5),'sample','FontSize',20)
for a=1:4
    ax(a).XTickLabel='';
    ax(a).FontSize=20;
end
ax(5).FontSize=20;
xlabel(ax(5),[xlabel_str '(seconds)'],'fontsize',20)    

%tscan   = 5
nb_seg=floor(timeaxis(end)/tscan);
if nb_seg<20
    error('too short. gather more data')
end


% sample rate channels
FS        = 1./nanmean(diff(timeaxis));
% number of samples per scan (1s) in channels
df        = 1/tscan;
f=(df:df:FS/2)'; % frequency vector for spectra
%% Length of the EPSI
T       = length(EPSI.epsitime);
%% define number of scan in the EPSI
Lscan   = tscan*2*f(end);
nbscan  = floor(T/Lscan);

%% we compute spectra on scan with 50% overlap
nbscan=2*nbscan-1;
channels=Meta_Data.PROCESS.channels;
nb_channels=length(channels);

%% split in segment of tscan length  the time series
indscan = arrayfun(@(x) (1+floor(Lscan/2)*(x-1):1+floor(Lscan/2)*(x-1)+Lscan-1),1:nbscan,'un',0);
nb_segment=numel(indscan);
%% select 2 segments around the middle of time serie
indscan1=indscan(floor(nb_segment/2)-10:floor(nb_segment/2)+10);
indexes=[indscan1{1}(1).*[1 1] indscan1{end}(end).*[1 1]];

%% plot on the time series the analysed block
for a=1:5
    hold(ax(a),'on')
    h(a)=fill(ax(a),timeaxis(indexes),[-2 2 2 -2],[.7 .7 .7]);
     h(a).FaceAlpha=.5;
end
legend(ax(1),{'a1','a2'})
legend(ax(2),{'a3'})
legend(ax(3),{'t1','t2'})
legend(ax(4),{'s1','s2'})
legend(ax(5),{'diff ramp'})
linkaxes(ax(1:5),'x');

%% clean time series
for c=1:nb_channels
    wh_channels=channels{c};
    EPSI.(wh_channels)=fillmissing(EPSI.(wh_channels)(indexes(1):indexes(3)),'linear');
    EPSI.(wh_channels)=filloutliers(EPSI.(wh_channels),'linear','movmedian',1000);
end

%% split data
data=zeros(nb_channels,numel(indscan1),Lscan);
for c=1:nb_channels
    wh_channels=channels{c};
    ind=find(cellfun(@(x) strcmp(x,wh_channels),channels));
    switch wh_channels
        case 't1'
            data(ind,:,:) = cell2mat(cellfun(@(x) EPSI.(wh_channels)(x-indexes(1)+1),indscan1,'un',0)).';
        case 't2'
            data(ind,:,:) = cell2mat(cellfun(@(x) EPSI.(wh_channels)(x-indexes(1)+1),indscan1,'un',0)).';
        case {'s1','s2'}
            data(ind,:,:) = cell2mat(cellfun(@(x) EPSI.(wh_channels)(x-indexes(1)+1),indscan1,'un',0)).';
        case {'a1','a2','a3'}
            data(ind,:,:) = cell2mat(cellfun(@(x) EPSI.(wh_channels)(x-indexes(1)+1),indscan1,'un',0)).';
    end
end

%% compute spectra
[f1,~,P11,~]=get_profile_spectrum(data,f);
indf1=find(f1>=0);
indf1=indf1(1:end-1);
f1=f1(indf1);
P11= 2*P11(:,:,indf1);

%% get TF
h_freq=get_filters_MADRE(Meta_Data,f1);

%% correct the TF
a1m=squeeze(nanmean(P11(6,:,:),2))./h_freq.electAccel(:);
a2m=squeeze(nanmean(P11(7,:,:),2))./h_freq.electAccel(:);
a3m=squeeze(nanmean(P11(8,:,:),2))./h_freq.electAccel(:);

s1m=squeeze(nanmean(P11(3,:,:),2))./h_freq.shear(:);
s2m=squeeze(nanmean(P11(4,:,:),2))./h_freq.shear(:);

t1m=squeeze(nanmean(P11(1,:,:),2))./h_freq.electFPO7(:);
t2m=squeeze(nanmean(P11(2,:,:),2))./h_freq.electFPO7(:);

%% get noise
Fn    = .5*FS;  % Nyquist frequency
FR    = 2.5;    % Full range in Volts
def_noise=@(x)((FR/2^x)^2 /Fn);
Accelnoise=45e-6^2+0*f1;

% polyfit the data
logf=log10(f1);
Empnoise=log10(squeeze(nanmean(P11(1,:,:),2)));
Empnoiseshear=log10(squeeze(nanmean(P11(3,:,:),2)));
Emp_FPO7noise=polyfit(logf(2:end).',Empnoise(2:end),3);
Emp_shearnoise=polyfit(logf(2:end).',Empnoiseshear(2:end),3);
%test_noise=polyval(Emp_FPO7noise,logf);
n3=Emp_FPO7noise(1);
n2=Emp_FPO7noise(2);
n1=Emp_FPO7noise(3);
n0=Emp_FPO7noise(4);

n3s=Emp_shearnoise(1);
n2s=Emp_shearnoise(2);
n1s=Emp_shearnoise(3);
n0s=Emp_shearnoise(4);

test_noise=n0+n1.*logf+n2.*logf.^2+n3.*logf.^3;
test_snoise=n0s+n1s.*logf+n2s.*logf.^2+n3s.*logf.^3;

%% plot spectra
ax(6)=subplot('Position',[.1 .05 .8 .52]);

hold(ax(6),'on')
l0=loglog(ax(6),f1,squeeze(nanmean(P11(1,:,:),2)),'--','Color',cmap(4,:));
loglog(ax(6),f1,squeeze(nanmean(P11(2,:,:),2)),'--','Color',cmap(5,:))
loglog(ax(6),f1,squeeze(nanmean(P11(3,:,:),2)),'--','Color',cmap(6,:))
loglog(ax(6),f1,squeeze(nanmean(P11(4,:,:),2)),'--','Color',cmap(7,:))
loglog(ax(6),f1,squeeze(nanmean(P11(6,:,:),2)),'--','Color',cmap(1,:))
loglog(ax(6),f1,squeeze(nanmean(P11(7,:,:),2)),'--','Color',cmap(2,:))
loglog(ax(6),f1,squeeze(nanmean(P11(8,:,:),2)),'--','Color',cmap(3,:))

% 
l1=loglog(ax(6),f1,t1m,'d-','Color',cmap(4,:));
l2=loglog(ax(6),f1,t2m,'p-','Color',cmap(5,:));
l3=loglog(ax(6),f1,s1m,'d-','Color',cmap(6,:));
l4=loglog(ax(6),f1,s2m,'p-','Color',cmap(7,:));
l5=loglog(ax(6),f1,a1m,'d-','Color',cmap(1,:));
l6=loglog(ax(6),f1,a2m,'p-','Color',cmap(2,:));
l7=loglog(ax(6),f1,a3m,'s-','Color',cmap(3,:));
set(ax(6),'Xscale','log','Yscale','log')


% bit noise
n20=loglog(ax(6),f1,f1*0+def_noise(20),'--','Color',[.5 .5 .5],'linewidth',2);
n24=loglog(ax(6),f1,f1*0+def_noise(24),'--','Color',[.1 .1 .1],'linewidth',2);
n16=loglog(ax(6),f1,f1*0+def_noise(16),'.-','Color',[.3 .3 .3],'linewidth',2);
An=loglog(ax(6),f1,Accelnoise,'--','Color',[.1 .1 .1],'linewidth',2);
Emp=loglog(ax(6),f1,10.^test_noise,'m-','linewidth',2);
Emps=loglog(ax(6),f1,10.^test_snoise,'c-','linewidth',2);

grid(ax(6),'on')
legend([l0,l1 l2 l3 l4 l5 l6 l7 n24 n20 n16 An Emp Emps],{'no TF','t1','t2','s1','s2','a1','a2','a3','24 bit','20 bit','16 bit','Accel noise','TF-noise','Sh-noise'},'location','SouthWest')
set(ax(6),'fontsize',30)
ylabel(ax(6),'V^2 / Hz','fontsize',30)
xlabel(ax(6),'Hz','fontsize',30)
ax(6).XLim=[1/tscan f(end)];
title(ax(1),[Meta_Data.MADRE.rev '-' Meta_Data.MADRE.SN '-' Meta_Data.MAP.rev '-' Meta_Data.MAP.SN],'fontsize',25)


fig=gcf;fig.PaperPosition = [0 0 30 30];
print(fullfile(Meta_Data.Epsipath,[Meta_Data.deployment '.png']),'-dpng')

Answer1=input('Do you want to save?(y/n)','s');
while 1
    switch Answer1
        case 'y'
            eval(['!mkdir ' Meta_Data.CALIpath])
             print(fullfile(Meta_Data.CALIpath,[Meta_Data.deployment '.png']),'-dpng')
             save(fullfile(Meta_Data.CALIpath,[Meta_Data.deployment '-FPO7_noise.mat']),'n0','n1','n2','n3')
             save(fullfile(Meta_Data.CALIpath,[Meta_Data.deployment '-shear_noise.mat']),'n0s','n1s','n2s','n3s')
             save(fullfile(Meta_Data.CALIpath,[Meta_Data.deployment '-spectrum_calib.mat']),'f1','P11')
             save(fullfile(Meta_Data.CALIpath, ...
            ['Meta_' Meta_Data.mission '_' Meta_Data.deployment '.mat']),'Meta_Data')
            break;
        case 'n'
            break;
        otherwise
            Answer1=input('Do you want to save?(y/n)','s');
    end
end

