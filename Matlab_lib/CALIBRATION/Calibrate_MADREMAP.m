%function Calibrate_MADREMAP(Meta_Data,file)

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




addpath ../EPSILON/toolbox/

file='/Users/aleboyer/ARNAUD/SCRIPPS/EPSILOMETER/CALIBRATION/ELECTRONICS/MADRE_MAP/MM1/epsi/MM1.mat'
%file='/Volumes/DataDrive/GRANITE/For_Arnaud_Sept2018/benchtest_for_alb/aug30_11/epsi/epsi_aug30_11.mat'
load('/Users/aleboyer/ARNAUD/SCRIPPS/EPSILOMETER/CALIBRATION/ELECTRONICS/MADRE_MAP/MM1/raw/Meta_DEV_MM1.mat')

EPSI=load(file);
F=fieldnames(EPSI);
T1=numel(EPSI.epsitime);
for f=1:length(F)
    wh_f=F{f};
    if numel(EPSI.(wh_f))==numel(EPSI.epsitime)
        EPSI.(wh_f)=EPSI.(wh_f)(T1/2:end);
    end
end

tscan     =  5;
% sample rate channels
FS        = 325;
% number of samples per scan (1s) in channels
df        = 1/tscan;
f=(df:df:FS/2)'; % frequency vector for spectra
%% Length of the EPSI
T       = length(EPSI.epsitime);
df      = f(1);
%% define number of scan in the EPSI
Lscan   = tscan*2*f(end);
nbscan  = floor(T/Lscan);

%% we compute spectra on scan with 50% overlap
nbscan=2*nbscan-1;
channels=strsplit('t1,t2,s1,s2,c,a1,a2,a3',',');
nb_channels=length(channels);


indscan = arrayfun(@(x) (1+floor(Lscan/2)*(x-1):1+floor(Lscan/2)*(x-1)+Lscan-1),1:nbscan,'un',0);

%TODO probably check nan at previous step
All_channels=fields(EPSI);
for c=1:length(All_channels)
    wh_channels=All_channels{c};
    EPSI.(wh_channels)=fillmissing(EPSI.(wh_channels),'linear');
    EPSI.(wh_channels)=filloutliers(EPSI.(wh_channels),'linear');
end


data=zeros(nb_channels,nbscan,Lscan);
for c=1:length(All_channels)
    wh_channels=All_channels{c};
    ind=find(cellfun(@(x) strcmp(x,wh_channels),channels));
    switch wh_channels
        case 't1'
            data(ind,:,:) = cell2mat(cellfun(@(x) EPSI.(wh_channels)(x),indscan,'un',0).');
        case 't2'
            data(ind,:,:) = cell2mat(cellfun(@(x) EPSI.(wh_channels)(x),indscan,'un',0).');
        case {'s1','s2'}
            data(ind,:,:) = cell2mat(cellfun(@(x) EPSI.(wh_channels)(x),indscan,'un',0).');
        case {'a1','a2','a3'}
            data(ind,:,:) = cell2mat(cellfun(@(x) EPSI.(wh_channels)(x),indscan,'un',0).');
    end
end


[f1,~,P11,~]=get_profile_spectrum(data,f);
indf1=find(f1>=0);
indf1=indf1(1:end-1);
f1=f1(indf1);
P11= 2*P11(:,:,indf1);


h_freq=get_filters_MADRE(Meta_Data,f1);


a1m=squeeze(nanmean(P11(6,:,:),2))./h_freq.electAccel(:);
a2m=squeeze(nanmean(P11(7,:,:),2))./h_freq.electAccel(:);
a3m=squeeze(nanmean(P11(8,:,:),2))./h_freq.electAccel(:);


s1m=squeeze(nanmean(P11(3,:,:),2))./h_freq.shear(:);
s2m=squeeze(nanmean(P11(4,:,:),2))./h_freq.shear(:);

t1m=squeeze(nanmean(P11(1,:,:),2))./h_freq.electFPO7(:);
t2m=squeeze(nanmean(P11(2,:,:),2))./h_freq.electFPO7(:);

Fn    = .5*FS;  % Nyquist frequency
FR    = 2.5;    % Full range in Volts
def_noise=@(x)((FR/2^x)^2 /Fn);
Accelnoise=45e-6^2+0*f1;


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


figure
hold on
l0=loglog(f1,squeeze(nanmean(P11(1,:,:),2)),'m--');
loglog(f1,squeeze(nanmean(P11(2,:,:),2)),'m--')
loglog(f1,squeeze(nanmean(P11(3,:,:),2)),'b--')
loglog(f1,squeeze(nanmean(P11(4,:,:),2)),'b--')
loglog(f1,squeeze(nanmean(P11(6,:,:),2)),'r--')
loglog(f1,squeeze(nanmean(P11(7,:,:),2)),'r--')
loglog(f1,squeeze(nanmean(P11(8,:,:),2)),'r--')
set(gca,'Xscale','log','Yscale','log')
hold on
l1=loglog(f1,t1m,'md-');
l2=loglog(f1,t2m,'mp-');
l3=loglog(f1,s1m,'bd-');
l4=loglog(f1,s2m,'bp-');
l5=loglog(f1,a1m,'rd-');
l6=loglog(f1,a2m,'rp-');
l7=loglog(f1,a3m,'rs-');
set(gca,'Xscale','log','Yscale','log')

n20=loglog(f1,f1*0+def_noise(20),'--','Color',[.5 .5 .5],'linewidth',2);
n24=loglog(f1,f1*0+def_noise(24),'--','Color',[.1 .1 .1],'linewidth',2);
n16=loglog(f1,f1*0+def_noise(16),'.-','Color',[.3 .3 .3],'linewidth',2);
An=loglog(f1,Accelnoise,'--','Color',[.1 .1 .1],'linewidth',2);
Emp=loglog(f1,10.^test_noise,'m-','linewidth',2);
Emps=loglog(f1,10.^test_snoise,'c-','linewidth',2);

grid on
legend([l0,l1 l2 l3 l4 l5 l6 l7 n24 n20 n16 An Emp Emps],{'no TF','t1','t2','s1','s2','a1','a2','a3','24 bit','20 bit','16 bit','Accel noise','TF-noise','Sh-noise'},'location','Southwest')
set(gca,'fontsize',15)
ylabel('V^2 / Hz','fontsize',15)
xlabel('Hz','fontsize',15)
title([ Meta_Data.MADRE.rev '-' Meta_Data.MADRE.SN '-' Meta_Data.MAP.rev '-' Meta_Data.MAP.SN],'fontsize',25)


fig=gcf;fig.PaperPosition = [0 0 20 10];
print([Meta_Data.L1path Meta_Data.deployement '.png'],'-dpng')

save([Meta_Data.CALIpath 'FPO7_noise.mat'],'n0','n1','n2','n3')
save([Meta_Data.CALIpath 'shear_noise.mat'],'n0s','n1s','n2s','n3s')

