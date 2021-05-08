fid=fopen('charge_amp_08012019.txt');
C=textscan(fid,'%f %f,%f %f,%f','Headerlines',1);

Vin=complex(C{4},C{5});
Vout=complex(C{2},C{3});

coef_filt=abs(Vout./Vin);
freq=C{1};

figure
hold on
semilogx(C{1},coef_filt)
set(gca,'Xscale','log')
pause
close all
% save('ca_01082019_TF.mat','coef_filt','freq');



