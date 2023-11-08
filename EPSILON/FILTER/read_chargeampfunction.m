fid=fopen('shearprobeAMP_v6.txt');
fid=fopen('shearprobeAMP_v6.txt');
C=textscan(fid,'%s %s','Headerlines',1);
cfreq=char(C{1});
for i=1:length(cfreq)
     freq(i)=str2double(cfreq(i,:));
end

value=C{2}(:);
value1=zeros(1,length(value));
for n=1:length(value)
    if value{n}(2)=='-'
        
        value1(n)=str2num(value{n}(2:23));
    else
        value1(n)=str2num(value{n}(2:22));
    end
end

coef_filt=10.^(value1/20); % value of the spice model given by sean are in dB
%save('charge_coeffilt.mat','coef_filt','freq');