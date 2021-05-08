function ax=MHA_diagnostic1fcn(MS,EpsiProfiles,CTDProfiles,Meta_Data,wh,z,N,nfft,nfftc,up_down)
%function ax=MHA_diagnostic1fcn(MS,EpsiProfiles,Meta_Data,wh,z,N,nfft,nfftc,up_down)
%% Plot a 4-panel plot of a) Arnaud's calculated shear spectra v freq, b) vel spectra v freq; c) accel spectra and d) coherence v freq.
% The last three are computed for profile wh with +/-N points centered around the
% specified bin (at depth z), with windows of length nfft for spectra and
% nfftc for coherence. Default values are nfft=4096,nfftc=1024, N=8000.
%
%A basic plot of fall rate, epsilon, T and S is also included.
%
%MHA Nov 2019

if nargin < 9
    up_down='down'; %default for down-data (epsi-fish)
end

if nargin < 8
    nfftc=1024;
end
if nargin < 7 
    nfft=4096;
end
if nargin < 6
    N=8000;
end

if ~isfield(Meta_Data,'SN')
    Meta_Data.SN=Meta_Data.vehicle_name;
end

samplefreq=325; %epsi sample freq
whsh=1:2;

%%
%Plot one spectrum... this is the shear wavenumber spectrum plotted versus
%freq

fig1=figure(1);
clf

xl=[1e-1 200];

ax=MySubplot(.1,.1,0.02,.1,.1,.02,2,2);
setpp(7,7)
axes(ax(1))
%whz=min(find(MS{wh}.pr > z));
whz=find(abs(MS{wh}.pr - z) ==min(abs(MS{wh}.pr -z)));

whsh=1:2;
h=loglog(MS{wh}.k.*MS{wh}.w(whz),squeeze(MS{wh}.Pshear_k(whz,:,whsh)),'r-',...
    MS{wh}.k.*MS{wh}.w(whz),squeeze(MS{wh}.Ppan(whz,:,whsh)),'r--');shg
% ALB lc might line color?
% lc(h([1 3]),'b')
h(1).Color='b';
h(3).Color='b';

xlim(xl)
legend(h,'v1','v2','location','southwest')
ylabel('\Phi_{shear}(\omega) / s^{-1}/cpm')
%xlabel('Hz')
grid
SubplotLetter('Computed spectra from MS structure',.01,.9)
title([Meta_Data.mission ', SN ' Meta_Data.SN ', dep ' Meta_Data.deployment ', profile #' num2str(wh) ', z=' num2str(MS{wh}.pr(whz)) ', \epsilon_1=' num2str(MS{wh}.epsilon(whz,1)) ', \epsilon_2=' num2str(MS{wh}.epsilon(whz,2))] )
xtloff

%%


% Now compute my own shear, accel spectra and coherences

%load('/Users/malford/GoogleDrive/Data/epsi/PISTON/Cruises/sr1914/data/epsi/PC2/d10/L1/Profiles_d10.mat')

%wh=6; %OK this is the range of data for a particular epsilon scan

 %Use the same data
iv=MS{wh}.indscan{whz};

iv=fix(mean(MS{wh}.indscan{whz}))+(-N:N); 

%Get the data - for vel and accel.  These should be in physical units.
%up_down='up';
whfield=['data' up_down];
v1=detrend(EpsiProfiles.(whfield){wh}.s1(iv));
v2=detrend(EpsiProfiles.(whfield){wh}.s2(iv));
a1=detrend(EpsiProfiles.(whfield){wh}.a1(iv));
if isfield(EpsiProfiles.(whfield){wh},'a2')
    a2=detrend(EpsiProfiles.(whfield){wh}.a2(iv)); %The other acceleration channel
else
    a2=nan*a1;
end    
a3=detrend(EpsiProfiles.(whfield){wh}.a3(iv)); %The other acceleration channel

%nfft=512;
[Pa1,fe] = pwelch(a1,nfft,[],nfft,samplefreq,'psd');
[Pa2,fe] = pwelch(a2,nfft,[],nfft,samplefreq,'psd');
[Pa3,fe] = pwelch(a3,nfft,[],nfft,samplefreq,'psd');
[Pv1,fe] = pwelch(v1,nfft,[],nfft,samplefreq,'psd'); 
[Pv2,fe] = pwelch(v2,nfft,[],nfft,samplefreq,'psd');


% Now compute k
ke=fe./MS{wh}.w(whz);

%% Plot the spectra
%figure(3)

%clf
%loglog(fe,Pe);shg
axes(ax(2))
loglog(fe,Pv1,'b-',fe,Pv2,'r');shg
grid
%freqline(50)
set(gca,'yaxislocation','right')
ylabel('s^{-2}/Hz')
title('Velocity frequency spectrum')
xtloff
xlim(xl)
SubplotLetter(['Velocity spectra, N=' num2str(N) ', nfft=' num2str(nfft)],.01,.9)

%% accel spectra
Accelnoise=45e-6^2+0*MS{wh}.f;
Accelnoise_mmp=35e-6^2+0*MS{wh}.f;
accel_chindex1=4:5; %6/7/8 but for NISKINE we were only transmitting 5 channels. 
accel_chindex2=5; %6/7/8 but for NISKINE we were only transmitting 5 channels. 
%figure(4)
%clf
%loglog(fe,Pe);shg
axes(ax(3))
loglog(fe,Pa1,'r-',fe,Pa2,'g-',fe,Pa3,'b-',MS{wh}.f,Accelnoise,'k--');shg;grid
xlim(xl)
title('Acceleration frequency spectrum')
ylabel('g^2/Hz')
legend('a1','a2','a3','noise spec','location','southwest')
xlabel('Hz')
SubplotLetter(['Acceleration spectra, N=' num2str(N) ', nfft=' num2str(nfft)],.01,.9)

%% Coherence
%[Cxy,F] = mscohere(X,Y,WINDOW,NOVERLAP,F,Fs)
%data1e=detrend(v1e);
%data2e=detrend(a1e);
[Cv1a1,fCe] = mscohere(v1,a1,nfftc,[],nfftc,samplefreq);
[Cv2a1,fCe] = mscohere(v2,a1,nfftc,[],nfftc,samplefreq);
[Cv1a2,fCe] = mscohere(v1,a2,nfftc,[],nfftc,samplefreq);
[Cv2a2,fCe] = mscohere(v2,a2,nfftc,[],nfftc,samplefreq);
[Cv1a3,fCe] = mscohere(v1,a3,nfftc,[],nfftc,samplefreq);
[Cv2a3,fCe] = mscohere(v2,a3,nfftc,[],nfftc,samplefreq);
%samplefreq=400;
%figure(5);clf;
axes(ax(4))
h=semilogx(fCe,Cv1a1,'b-',fCe,Cv2a1,'r-',fCe,Cv1a2,'b--',fCe,Cv2a2,'r--',fCe,Cv1a3,'b-',fCe,Cv2a3,'r')
% 
% lw(h(5:6),2)
h(5).LineWidth=2;
h(6).LineWidth=2;

xlim(xl)
grid
title('Coherence')
SubplotLetter(['Coherence, N=' num2str(N) ', nfft=' num2str(nfftc)],.01,.9)

legend('v1a1','v2a1','v1a2','v2a2','v1a3','v2a3')
set(gca,'yaxislocation','right')
ylabel('Coherence')
xlabel('Hz')

%%
fig2=figure(2);

clf
ax=MySubplot(.1,.1,0,.1,.1,0,4,1);
axes(ax(1))
plot(MS{wh}.w,MS{wh}.pr);axis ij; shg
xlabel('w / m~s^{-1}')
grid
axes(ax(2))
%semilogx(MS{wh}.epsilon(:,whsh),MS{wh}.pr,MS{wh}.epsilon_co(:,whsh),MS{wh}.pr);axis ij; shg
semilogx(MS{wh}.epsilon(:,whsh),MS{wh}.pr);axis ij; shg
% ALB what is ytloff
% ytloff
xlim([1e-11 1e-6])
grid
hold on
tmp=xlim;
plot(tmp,[z z],'k--')
plot(tmp,[z z]+N./samplefreq.*MS{wh}.w(whz),'r--')
plot(tmp,[z z]-N./samplefreq.*MS{wh}.w(whz),'r--')

hold off
xlabel('\epsilon / W~kg^{-1}')
axes(ax(3))
plot(MS{wh}.t,MS{wh}.pr);axis ij; shg
% % ALB what is ytloff
% ytloff
xlabel('T / ^oC')
grid
axes(ax(4))
plot(MS{wh}.s,MS{wh}.pr);axis ij; shg
% % ALB what is ytloff
% ytloff
xlabel('S/ psu')
grid
linkaxes(ax,'y')

%Get back to the spectra plot
figure(1)
fig1.PaperPosition=[0 0 18 12];
figname=[Meta_Data.mission '_SN' Meta_Data.SN '_dep_' Meta_Data.deployment '_profile' num2str(wh) '_z' num2str(MS{wh}.pr(whz)) '_spectra.png'];
print('-dpng2',figname)

figure(2)
fig2.PaperPosition=[0 0 15 15];

figname=[Meta_Data.mission '_SN' Meta_Data.SN '_dep_' Meta_Data.deployment '_profile' num2str(wh) '_z' num2str(MS{wh}.pr(whz)) '_WTS.png'];
print('-dpng2',figname)

figure(1)






