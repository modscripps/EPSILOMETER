fid=fopen('tcurve_shr_chrgamp_waa.txt');
C=textscan(fid,'%f %f %f %f %f','Headerlines',1);

freq=C{1};

value=C{2};
coef_filt=10.^(value/20); % value of the spice model given by sean are in dB
coef_filt1=10.^(value/10); % value of the spice model given by sean are in dB

semilogx(freq,coef_filt.^2,'r')
hold on
semilogx(freq,coef_filt1.^2,'b')

save('shr_chramp_waa.mat','coef_filt','freq');