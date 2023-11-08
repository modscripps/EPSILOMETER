% tdiff pre tdiff mmp noise
mmpTdiff=load ('~/ARNAUD/SCRIPPS/EPSILOMETER/EPSILON/toolbox/FILTER/mmp_tdiff.mat','f','helectronics');
Tdiff_filter = load('FILTER/Tdiff_filt.mat','coef_filt','freq');


FPO7noise=load('~/ARNAUD/SCRIPPS/EPSILOMETER/CALIBRATION/ELECTRONICS/FPO7_noise.mat','n0','n1','n2','n3');

Meta_Data.epsi.s1.ADCfilter='sinc4';
Meta_Data.epsi.t1.ADCfilter='sinc4';
Meta_Data.epsi.a1.ADCfilter='sinc4';
Meta_Data.MAP.temperature='Tdiff';
f=1/3:1/3:320;
H=get_filters_MADRE(Meta_Data,f);
epsiTF=H.electFPO7.^2 .* H.Tdiff.^2;

logf=log10(f);
epsi_noise=10.^(FPO7noise.n0+FPO7noise.n1.*logf+FPO7noise.n2.*logf.^2+FPO7noise.n3.*logf.^3);


mmp_n0=-7.8; mmp_n1=-0.0634538; mmp_n2=0.3421899; mmp_n3=-0.3141283;
logf=log10(mmpTdiff.f);
mmp_noise=10.^(mmp_n0+mmp_n1.*logf+mmp_n2.*logf.^2+mmp_n3.*logf.^3);

% beta in mmp
alpha=-0.05;
R1=1e4;
R2=1e6;
Rt=820e3;
E0=10;
%beta_mmp=3.7097e-29;
beta=alpha*(R1+R2)*Rt*E0./(R1+R2+Rt).^2;
Gtl=2.529;

kb = 1.3806e-23; % boltzman constant
T10  = 10+273;  % 10 degres C
R750  = 750e3;   % 150 kohms
R150  = 150e3;   % 150 kohms
Df = 1;
johnson150=4*kb*T10*R150*Df./(1/3);
johnson750=4*kb*T10*R750*Df./(1/3);


%%
close all

F1=figure;
loglog(mmpTdiff.f,mmp_noise,'b-.','linewidth',2)
hold on
loglog(mmpTdiff.f,sqrt(mmpTdiff.helectronics),'b--','linewidth',2)
loglog(mmpTdiff.f,mmp_noise./mmpTdiff.helectronics/beta.^2/Gtl.^2,'b','linewidth',3)


loglog(f,epsi_noise,'r-.','linewidth',2)
loglog(f,sqrt(epsiTF),'r--','linewidth',2)
loglog(f,epsi_noise./epsiTF*30^2,'r','linewidth',3)
loglog(f,f*0+johnson750,'Color',[.5 .5 .5],'linewidth',2)
loglog(f,f*0+johnson150,'Color',[.2 .2 .2],'linewidth',2)
legend('MMP noise','MMP TF (no units)','MMP noise./TF','EPSI noise', ...
       'EPSI TF (no units)','EPSI noise/TF', ...
       'Johnson 10 C 750 kOhms','Johnson 10 C 150 kOhms','location','southwest')

setaxis_frequency_spectra('Hz','T^2/Hz');
fig=gcf;fig.PaperPosition=[0 0 20 20];
print('MMP_EPSI_noise_comparison.png','-dpng2')

