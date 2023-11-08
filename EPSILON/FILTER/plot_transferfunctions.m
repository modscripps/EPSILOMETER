nb_headerline=6;
rootpath='/Volumes/MOD_dev/Projects/Epsilometer/CALIBRATION/AFE/';
channels={'t1','t2','s1','s2'};
for i=[2 3 4 6 8]
    for c=1:length(channels)
        wh_ch=channels{c};
        filestr=sprintf('SN%i/MAP_SN%i',i,i);
        SNstr=sprintf('SN%i',i);
        filename=sprintf('%s%s_%s.csv',rootpath,filestr,wh_ch);
        TF.(SNstr).(wh_ch)=read_digilient_network_analysis(filename);
    end
end

%% current trnasfet functions

fid=fopen('/Users/aleboyer/ARNAUD/SCRIPPS/EPSILOMETER/EPSILON/toolbox/FILTER/EpsidTdt_H/Epsi_FPO7dTdt_ri.txt');
C=textscan(fid,'%f %f,%f %f,%f','Headerlines',1);
Vin=complex(C{2},C{3});
Vout=complex(C{4},C{5});

Temp.coef_filt=abs(Vout./Vin);
Temp.freq=C{1};


Charge_Amp=load('charge_coeffilt.mat');


%% plots


close all; 
figure(1)
cmap=colormap(jet(5));
hold on
count=0
ll=[];
for i=[2 3 4 6 8]
    count=count+1;
    SNstr=sprintf('SN%i',i);
    legendstr{count}=SNstr;
    l=semilogx(TF.(SNstr).t1.f,TF.(SNstr).t1.Vout./TF.(SNstr).t1.Vin,'Color',cmap(count,:),'linewidth',1.5);
    semilogx(TF.(SNstr).t2.f,TF.(SNstr).t2.Vout./TF.(SNstr).t2.Vin,'Color',cmap(count,:),'linewidth',1.5)
    ll=[ll l];
end
l=semilogx(Temp.freq,Temp.coef_filt,'k')
ll=[ll l];
legendstr{count+1}='SPICE';
title('Temp','fontsize',20)
legend(ll,legendstr)
grid on
xlim([1e-2 1e4])
set(gca,'fontsize',20)
set(gca,'Xscale','log')
ylabel('H (no units)','fontsize',20)
xlabel('Hz','fontsize',20)
fig=gcf;fig.PaperPosition=[0 0 15 10];
print('-dpng2','~/ARNAUD/SCRIPPS/EPSILOMETER/CALIBRATION/ELECTRONICS/AFE_TF_temperature.png')


%% shear


close all; 
figure(1)
cmap=colormap(jet(5));
hold on
count=0
ll=[];
for i=[2 3 4 6 8]
    count=count+1;
    SNstr=sprintf('SN%i',i);
    legendstr{count}=SNstr;
    l=semilogx(TF.(SNstr).t1.f,20*log10(TF.(SNstr).s1.Vout./TF.(SNstr).s1.Vin),'Color',cmap(count,:),'linewidth',1.5);
    semilogx(TF.(SNstr).t2.f,20*log10(TF.(SNstr).s2.Vout./TF.(SNstr).s2.Vin),'Color',cmap(count,:),'linewidth',1.5)
    ll=[ll l];
end
l=semilogx(Charge_Amp.freq,20*log10(Charge_Amp.coef_filt),'b')

ll=[ll l];
legendstr{count+1}='SPICE';
title('Shear','fontsize',20)
legend(ll,legendstr)
%ylim([0 1])
grid on
xlim([1e-2 1e4])
set(gca,'fontsize',20)
set(gca,'Xscale','log')
%ylabel('H (no units)','fontsize',20)
ylabel('dB  (no units)','fontsize',20)
xlabel('Hz','fontsize',20)
fig=gcf;fig.PaperPosition=[0 0 15 10];
print('-dpng2','~/ARNAUD/SCRIPPS/EPSILOMETER/CALIBRATION/ELECTRONICS/AFE_TF_shear_dB.png')

% freq=Charge_Amp.freq;
% %save('charge_coeffilt_09312019.mat','freq','coef_filt');




