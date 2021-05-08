function [F1,F2,legend_string]=plot_binned_chi(Chi_class,Meta_Data,title_string,indplot)
%  input: 
%  Chi_class : averaged temperature gradiant spectra classified by chi value
%  Meta_Data : path to calibration file and EPSI configuration needed to process the
%              epsi data
%  title_string: title of the plot
%  indplot: index of the chi value you want to plot (TODO: directly the chi values instead of the indexes)
%
%  Created by Arnaud Le Boyer on 7/28/18.



F1 = figure(1);clf
set(F1,'position',[100 50 2000 1100])
axes('position',[.05 .1 .82 .78])

L=size(Chi_class.k,2);
[Nchi,Neps]=size(Chi_class.nbin11);

if nargin<4
    indplot=[1 Nchi*Neps];
end

NOISE=load('toolbox/PLOTS/comparison_temp_granite_sproul.mat');
H=get_filters_MADRE(Meta_Data,NOISE.k_granite);
noise=NOISE.spec_granite./H.FPO7(.6);
noise=interp1(NOISE.k_granite/.6,noise,Chi_class.k);
% 
% dataP=reshape(Chi_class.mPphiT11,[Nchi*Neps,L]);
% inddata=find(nansum(dataP,2)>0);
% indplot1=max([indplot(1) 1]):min([indplot(end) length(inddata)]);
% inddata=inddata(indplot1);
% [chi,epsilon]=meshgrid(Chi_class.chibin,Chi_class.epsibin);
% chi=chi(:);epsilon=epsilon(:);
% databatch=reshape(Chi_class.Pbatch11,[Nchi*Neps,L]);
% databatch(databatch==0)=nan;
% cmap = parula(length(inddata)+1);
% 
% bb = loglog(Chi_class.kbatch,databatch(inddata,:),'color',[1 1 1]*.7,'linewidth',1);
% hold on
% b1 = loglog(Chi_class.kbatch,squeeze(Chi_class.Pbatch11(1,1,:)),'o-','color',[1 1 1]*.25,'linewidth',2);
% b2 = loglog(Chi_class.kbatch,squeeze(Chi_class.Pbatch11(1,end,:)),'+-','color',[1 1 1]*.4,'linewidth',2);
% b3 = loglog(Chi_class.kbatch,squeeze(Chi_class.Pbatch11(end,1,:)),'d-','color',[1 1 1]*.5,'linewidth',2);
% b4 = loglog(Chi_class.kbatch,squeeze(Chi_class.Pbatch11(end,end,:)),'s-','color',[1 1 1]*.1,'linewidth',2);
% for i=1:length(inddata)
%     j=inddata(i);
%     ll(i)=loglog(Chi_class.k,squeeze(dataP(j,:)),'Color',cmap(i,:),'linewidth',2);
% end
% 
% lnoise=loglog(Chi_class.k,(2*pi*Chi_class.k).^2 .* 62^2.*noise/.6,'r');
% 
% mink=Chi_class.k(find(~isnan(nansum(dataP,1)),1,'first'));
% maxk=Chi_class.k(find(~isnan(nansum(dataP,1)),1,'last'));
% minP=min(dataP(:));
% maxP=max(dataP(:));
% 
% 
% xlim([mink maxk])
% ylim([minP maxP])
% grid on
% xlabel('k [cpm]')
% ylabel('[$(\phi^T_k)^2$  / $(^{\circ}C . m^{-1})^2$ / cpm]','interpreter','latex')
% set(gca,'fontsize',20)
% title([title_string '- $\phi^{T22}_k$'],'interpreter','latex')
% 
% legend_string=arrayfun(@(x,y) sprintf('log_{10}(\\epsilon,\\phi)~%2.1f,%2.1f',x,y),...
%               log10(epsilon(inddata)),log10(chi(inddata)),'un',0);
% legend_string{length(inddata)+1}=sprintf('batchelor, log_{10}(\\epsilon,\\phi)~%2.1f,%2.1f',min(log10(chi)),min(log10(epsilon)));
% legend_string{length(inddata)+2}=sprintf('batchelor, log_{10}(\\epsilon,\\phi)~%2.1f,%2.1f',min(log10(chi)),max(log10(epsilon)));
% legend_string{length(inddata)+3}=sprintf('batchelor, log_{10}(\\epsilon,\\phi)~%2.1f,%2.1f',max(log10(chi)),min(log10(epsilon)));
% legend_string{length(inddata)+4}=sprintf('batchelor, log_{10}(\\epsilon,\\phi)~%2.1f,%2.1f',max(log10(chi)),max(log10(epsilon)));
% legend_string{length(inddata)+5}='current noise';
% 
% %hl=gridLegend([ll b1 b2 b3 b4],1,legend_string,'location','eastoutside');
% hl=legend([ll b1 b2 b3 b4 lnoise],legend_string,'location','eastoutside');
% set(hl,'fontsize',10)
% % 
% % axes('position',[.7 .1 .25 .82])
% % 
% % plot(log10(Chi_class.bin),Chi_class.nbin1,'-+k')
% % set(gca,'fontsize',20)
% % xlabel('log_{10} \epsilon')
% % xlim([min(log10(Chi_class.bin)) max(log10(Chi_class.bin))])
% % title('PDF \epsilon_1')
% 


F2 = figure(2);clf
clear ll
clear bb

set(F2,'position',[100 50 2000 1100])
axes('position',[.05 .1 .82 .78])

[Nchi,Neps,L]=size(Chi_class.Pbatch21);
dataP=reshape(Chi_class.mPphiT21,[Nchi*Neps,L]);
inddata=find(nansum(dataP,2)>0);
indplot2=max([indplot(1) 1]):min([indplot(end) length(inddata)]);
inddata=inddata(indplot2);
[chi,epsilon]=meshgrid(Chi_class.chibin,Chi_class.epsibin);
chi=chi(:);epsilon=epsilon(:);
databatch=reshape(Chi_class.Pbatch21,[Nchi*Neps,L]);
databatch(databatch==0)=nan;
cmap = parula(length(inddata)+1);

bb = loglog(Chi_class.kbatch,databatch(inddata,:),'color',[1 1 1]*.7,'linewidth',1);
hold on
b1 = loglog(Chi_class.kbatch,squeeze(Chi_class.Pbatch21(1,1,:)),'o-','color',[1 1 1]*.25,'linewidth',2);
b2 = loglog(Chi_class.kbatch,squeeze(Chi_class.Pbatch21(1,end,:)),'+-','color',[1 1 1]*.4,'linewidth',2);
b3 = loglog(Chi_class.kbatch,squeeze(Chi_class.Pbatch21(end,1,:)),'d-','color',[1 1 1]*.5,'linewidth',2);
b4 = loglog(Chi_class.kbatch,squeeze(Chi_class.Pbatch21(end,end,:)),'s-','color',[1 1 1]*.1,'linewidth',2);
for i=1:length(inddata)
    j=inddata(i);
    ll(i)=loglog(Chi_class.k,squeeze(dataP(j,:)),'Color',cmap(i,:),'linewidth',2);
end

lnoise=loglog(Chi_class.k,(2*pi*Chi_class.k).^2 .* 62^2.*noise/.6,'r');

mink=Chi_class.k(find(~isnan(nansum(dataP,1)),1,'first'));
maxk=Chi_class.k(find(~isnan(nansum(dataP,1)),1,'last'));
minP=min(dataP(:));
maxP=max(dataP(:));


xlim([mink maxk])
ylim([minP maxP])
grid on
xlabel('k [cpm]')
ylabel('[$(\phi^T_k)^2$  / $(^{\circ}C . m^{-1})^2$ / cpm]','interpreter','latex')
set(gca,'fontsize',20)
title([title_string '- $\phi^{T21}_k$'],'interpreter','latex')

legend_string=arrayfun(@(x,y) sprintf('log_{10}(\\epsilon,\\phi)~%2.1f,%2.1f',x,y),...
              log10(epsilon(inddata)),log10(chi(inddata)),'un',0);
legend_string{length(inddata)+1}=sprintf('batchelor, log_{10}(\\epsilon,\\phi)~%2.1f,%2.1f',min(log10(chi)),min(log10(epsilon)));
legend_string{length(inddata)+2}=sprintf('batchelor, log_{10}(\\epsilon,\\phi)~%2.1f,%2.1f',min(log10(chi)),max(log10(epsilon)));
legend_string{length(inddata)+3}=sprintf('batchelor, log_{10}(\\epsilon,\\phi)~%2.1f,%2.1f',max(log10(chi)),min(log10(epsilon)));
legend_string{length(inddata)+4}=sprintf('batchelor, log_{10}(\\epsilon,\\phi)~%2.1f,%2.1f',max(log10(chi)),max(log10(epsilon)));
legend_string{length(inddata)+5}='current noise';

%hl=gridLegend([ll b1 b2 b3 b4],1,legend_string,'location','eastoutside');
hl=legend([ll b1 b2 b3 b4 lnoise],legend_string,'location','eastoutside');
set(hl,'fontsize',10)






