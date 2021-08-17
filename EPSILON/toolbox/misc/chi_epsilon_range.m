ktemp=1.5e-7;
kvis=1e-6;
chilow=1e-12;
epsilonlow=1e-10;

[klow,Plow]=batchelor(epsilonlow,chilow,kvis,ktemp);

chihigh=3e-6;
epsilonhigh=1e-5;
[khigh,Phigh]=batchelor(epsilonhigh,chihigh,kvis,ktemp);

speed=.6;
M_Plow=max(Plow);
dTdtL=sqrt(M_Plow)*nanmean(diff(klow))./speed;
peak_kL=klow(Plow==M_Plow);
peak_fL=klow(Plow==M_Plow)*speed;
M_Phigh=max(Phigh);
dTdtH=sqrt(M_Phigh)*nanmean(diff(khigh))./speed;
peak_kH=khigh(Phigh==M_Phigh);
peak_fH=khigh(Phigh==M_Phigh)*speed;

%%
close all
loglog(klow,Plow,'b',khigh,Phigh,'r','linewidth',3);
grid on
hold on
blabla1={sprintf('kvis=%1.1e m^2 s^{-1} ',kvis),sprintf('ktemp=%1.1e m^2 s^{-1}',ktemp), ...
    sprintf('epsi=%1.1e W kg^{-1}',epsilonlow),sprintf('chi=%1.1e K^2 s^{-1}',chilow)};
blabla2={sprintf('kvis=%1.1e m^2 s^{-1}',kvis),sprintf('ktemp=%1.1e m^2 s^{-1}',ktemp), ...
    sprintf('epsi=%1.1e W kg^{-1}',epsilonhigh),sprintf('chi=%1.1e K^2 s^{-1}',chihigh)};
text(1,2e-7,blabla1,'Color','b','fontsize',20)
text(10,1e-4,blabla2,'Color','r','fontsize',20)

plot(peak_kH.*[1 1],[1e-10 M_Phigh],'r--','linewidth',3)
plot([1e-1 peak_kH],M_Phigh.*[1 1],'r--','linewidth',3)
plot(peak_kL.*[1 1],[1e-10 M_Plow],'b--','linewidth',3)
plot([1e-1 peak_kL],M_Plow.*[1 1],'b--','linewidth',3)

blabla3={sprintf('%3.1f cpm',peak_kH), ...
         sprintf('* %1.1f (m s^{-1}, fall rate)',speed),...
         sprintf('=%3.1f Hz',peak_fH)};
blabla4={sprintf('%3.1f cpm',peak_kL), ...
         sprintf('* %1.1f (m s^{-1}, fall rate)',speed),...
         sprintf('=%3.1f Hz',peak_fL)};
     
text(1.1*peak_kH,10^(mean(log10([1e-10 M_Phigh]))),blabla3,'Color','r','fontsize',20)
text(1.1*peak_kL,10^(mean(log10([1e-10 M_Plow]))),blabla4,'Color','b','fontsize',20)


blabla5={sprintf('%1.2e C^2 m^{-2} / cpm',M_Phigh),...
         sprintf('max dTdt= %2.2e C s^{-1}',dTdtH)};
blabla6={sprintf('%1.2e C^2 m^{-2} / cpm',M_Plow),...
         sprintf('max dTdt= %2.2e C s^{-1}',dTdtL)};

     
     
text(.2,.5*M_Phigh,blabla5,'Color','r','fontsize',20)
text(.2,.5*M_Plow,blabla6,'Color','b','fontsize',20)

xlim([1e-1 500])
ylim([1e-9 1e-0])
set(gca,'fontsize',20)  

xlabel('cpm','fontsize',20)
ylabel('C^2 m^{-2} / cpm','fontsize',20)
title('Batchelor','fontsize',20)

Tdiff_filter = load('~/ARNAUD/SCRIPPS/EPSILOMETER/EPSILON/toolbox/FILTER/Tdiff_filt.mat');
Tdiff_H = interp1(Tdiff_filter.freq,Tdiff_filter.coef_filt ,speed.*klow);
loglog(klow,1e-5*Tdiff_H,'--k')

blabla7={'10^{-5}*Mike''s Tdiff (no units)'};
text(1,1e-5,blabla7,'Color','k','fontsize',15)



%%
fig=gcf;fig.PaperPosition=[0 0 20 20];
print('-dpng2','/Users/aleboyer/ARNAUD/SCRIPPS/EPSILOMETER/EPSILON/toolbox/misc/chi_epsilon_range.png');

