function [Cu1a,Cu2a,sumCu1a,sumCu2a,fe]=mod_efe_scan_multivariate_coherence(scan,acceleration_channel,Meta_Data)
% Compute the coherence over the whole profile using all 3 axis (we could add more) 
% over the 1./tsan:Fs frequency frequency axis with nfft samples.

% written by aleboyer@ucsd.edu 10/01/2021


nfft=Meta_Data.PROCESS.nfft;
nfftc=Meta_Data.PROCESS.nfftc;
try
Fs=Meta_Data.PROCESS.Fs_epsi;
catch
   Fs=Meta_Data.AFE.FS; 
end
fc1=Meta_Data.PROCESS.fc1;
fc2=Meta_Data.PROCESS.fc2;

if isfinite(scan.s1_volt)
    [Cu1a,fe] = mscohere(detrend(scan.s1_volt),detrend(scan.(acceleration_channel)),nfftc,[],nfft,Fs);
    sumCu1a=sum(Cu1a(fe>fc1 & fe<fc2))*nanmean(diff(fe));
else
    fe = nan(nfft/2 + 1,1);
    Cu1a = nan(nfft/2 + 1,1);
    sumCu1a = nan;
end

if isfinite(scan.s2_volt)
    [Cu2a,~] = mscohere(detrend(scan.s2_volt),detrend(scan.(acceleration_channel)),nfftc,[],nfft,Fs);
    sumCu2a=sum(Cu2a(fe>fc1 & fe<fc2))*nanmean(diff(fe));
else
    Cu2a = nan(nfft/2 + 1,1);
    sumCu2a = nan;
end