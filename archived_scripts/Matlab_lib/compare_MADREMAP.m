function compare_MADREMAP(Meta_Data1,Meta_Data2)

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
P1=load(fullfile(Meta_Data1.CALIpath,[Meta_Data1.deployment '-spectrum_calib.mat']),'f1','P11');
P2=load(fullfile(Meta_Data2.CALIpath,[Meta_Data2.deployment '-spectrum_calib.mat']),'f1','P11');


%% get noise
Fn    = P1.f1(end);  % Nyquist frequency
FR    = 2.5;    % Full range in Volts
def_noise=@(x)((FR/2^x)^2 /Fn);
Accelnoise=45e-6^2+0*P1.f1;


%% plot spectra
close all
figure('units','inch','position',[0,0,35,15]);
cmap1=colormap(gray(15));
cmap2=colormap(parula(8));

ax(1)=axes('Position',[.1 .1 .8 .8]);

hold(ax(1),'on')
l11=loglog(ax(1),P1.f1,squeeze(nanmean(P1.P11(1,:,:),2)),'Color',[cmap1(6,:) .5],'linewidth',3);
l12=loglog(ax(1),P1.f1,squeeze(nanmean(P1.P11(2,:,:),2)),'Color',[cmap1(7,:) .5],'linewidth',3);
l13=loglog(ax(1),P1.f1,squeeze(nanmean(P1.P11(3,:,:),2)),'Color',[cmap1(4,:) .5],'linewidth',3);
l14=loglog(ax(1),P1.f1,squeeze(nanmean(P1.P11(4,:,:),2)),'Color',[cmap1(5,:) .5],'linewidth',3);
l15=loglog(ax(1),P1.f1,squeeze(nanmean(P1.P11(6,:,:),2)),'Color',[cmap1(1,:) .5],'linewidth',3);
l16=loglog(ax(1),P1.f1,squeeze(nanmean(P1.P11(7,:,:),2)),'Color',[cmap1(2,:) .5],'linewidth',3);
l17=loglog(ax(1),P1.f1,squeeze(nanmean(P1.P11(8,:,:),2)),'Color',[cmap1(3,:) .5],'linewidth',3);

l21=loglog(ax(1),P2.f1,squeeze(nanmean(P2.P11(1,:,:),2)),'Color',cmap2(6,:));
l22=loglog(ax(1),P2.f1,squeeze(nanmean(P2.P11(2,:,:),2)),'Color',cmap2(7,:));
l23=loglog(ax(1),P2.f1,squeeze(nanmean(P2.P11(3,:,:),2)),'Color',cmap2(4,:));
l24=loglog(ax(1),P2.f1,squeeze(nanmean(P2.P11(4,:,:),2)),'Color',cmap2(5,:));
l25=loglog(ax(1),P2.f1,squeeze(nanmean(P2.P11(6,:,:),2)),'Color',cmap2(1,:));
l26=loglog(ax(1),P2.f1,squeeze(nanmean(P2.P11(7,:,:),2)),'Color',cmap2(2,:));
l27=loglog(ax(1),P2.f1,squeeze(nanmean(P2.P11(8,:,:),2)),'Color',cmap2(3,:));

ax(1).XLim=[1/10 P1.f1(end)];
set(ax(1),'Xscale','log','Yscale','log')

% bit noise
n20=loglog(ax(1),P1.f1,P1.f1*0+def_noise(20),'--','Color',[.5 .5 .5],'linewidth',2);
n24=loglog(ax(1),P1.f1,P1.f1*0+def_noise(24),'--','Color',[.1 .1 .1],'linewidth',2);
n16=loglog(ax(1),P1.f1,P1.f1*0+def_noise(16),'.-','Color',[.3 .3 .3],'linewidth',2);
An=loglog(ax(1),P1.f1,Accelnoise,'--','Color',[.1 .1 .1],'linewidth',2);

grid(ax(1),'on')
legend([l11,l12 l13 l14 l15 l16 l17 ...
        l21,l22 l23 l24 l25 l26 l27 ...
        n24 n20 n16 An], ...
        {['t1-' Meta_Data1.deployment],'t2','s1','s2','a1','a2','a3', ...
         ['t1-' Meta_Data2.deployment],'t2','s1','s2','a1','a2','a3', ...
        '24 bit','20 bit','16 bit','Accel noise'},'location','SouthEastOutSide')
set(ax(1),'fontsize',30)
ylabel(ax(1),'V^2 / Hz','fontsize',30)
xlabel(ax(1),'Hz','fontsize',30)
title(ax(1),{[Meta_Data1.MADRE.rev '-' Meta_Data1.MADRE.SN '-' Meta_Data1.MAP.rev '-' Meta_Data1.MAP.SN], ...
             [Meta_Data2.MADRE.rev '-' Meta_Data2.MADRE.SN '-' Meta_Data2.MAP.rev '-' Meta_Data2.MAP.SN]},'fontsize',25)

fig=gcf;fig.PaperPosition = [0 0 30 20];
print(fullfile(Meta_Data2.CALIpath,[Meta_Data1.deployment '-' Meta_Data1.deployment '.png']),'-dpng')


% fig=gcf;fig.PaperPosition = [0 0 30 30];
% print(fullfile(Meta_Data.paths.profiles,[Meta_Data.deployment '.png']),'-dpng')

