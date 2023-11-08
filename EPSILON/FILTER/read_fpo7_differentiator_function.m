fid=fopen('EpsidTdt_H/Epsi_FPO7dTdt_ri.txt');
C=textscan(fid,'%f %f,%f %f,%f','Headerlines',1);

Vin=complex(C{2},C{3});
Vout=complex(C{4},C{5});

coef_filt=abs(Vout./Vin);
freq=C{1};

figure
hold on
semilogx(C{1},coef_filt)
set(gca,'Xscale','log')
pause
close all
save('Tdiff_filt.mat','coef_filt','freq');





% fid=fopen('EpsidTdt_H/Epsi_FPO7dTdt.txt');
% C=textscan(fid,'%s %s','Headerlines',1);
% cfreq=char(C{1});
% freq=zeros(length(cfreq),1);
% for i=1:length(cfreq)
%     freq(i)=str2double(cfreq(i,:));
% end
% 
% value=C{2}(:);
% value1=zeros(1,length(value));
% for n=1:length(value)
%     if value{n}(2)=='-'
%         
%         value1(n)=str2num(value{n}(2:23));
%     else
%         value1(n)=str2num(value{n}(2:22));
%     end
% end
% 
% coef_filt=10.^(value1/20); % value of the spice model given by sean are in dB

% save('Tdiff_filt.mat','coef_filt','freq');

