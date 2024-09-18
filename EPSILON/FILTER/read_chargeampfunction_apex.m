fid=fopen('cap1nFres200Meg_5KohmInput.txt');
C=textscan(fid,'%s %s %s','Headerlines',1);
fclose(fid);
cfreq=char(C{1});
for i=1:length(cfreq)
     freq(i)=str2double(cfreq(i,:));
end

db=cellfun(@str2num,C{2}(:));

semilogx(freq,db)
hold on
grid
plot(freq,-3+0.*freq,'k')

figure
coef_filt=10.^(db/20); % value of the spice model given by sean are in dB
semilogx(freq,coef_filt)


save('cap1nFres200Meg_5KohmInput.mat','coef_filt','freq');