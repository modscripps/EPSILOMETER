function [Cua,fe]=mod_efe_scan_shear_accel_coherence(scan,shear_channel,acceleration_channel,Meta_Data)
% get epsilon and the cutting frequency 


nfftc=Meta_Data.PROCESS.nfftc;
df=Meta_Data.df_epsi;


[Cua,fe] = mscohere(detrend(scan.(shear_channel)),detrend(scan.(acceleration_channel)),nfftc,[],nfftc,df);
