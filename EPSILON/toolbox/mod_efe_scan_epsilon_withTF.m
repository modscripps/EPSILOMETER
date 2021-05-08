function [P,Pv,Pvk,Psk,epsilon,fc,fe]=mod_efe_scan_epsilon_withTF(scan,vibration,shear_channel,Meta_Data)
% get epsilon and the cutting frequency 


nfft=Meta_Data.PROCESS.nfft;
Fs=Meta_Data.PROCESS.Fs_epsi;
fpump=Meta_Data.PROCESS.ctd_fc;

[P,fe] = pwelch(detrend(scan.(shear_channel)),nfft,[],nfft,Fs,'psd');

% get the filter transfer functions.
h_freq=Meta_Data.PROCESS.h_freq;

% remove the modeled signal created by vibration
Pp=P-vibration;
Pp(Pp<=0)=nan;
Pp=fillmissing(Pp,'linear');

filter_TF=(h_freq.shear .* haf_oakey(fe,scan.w));
Pv   = Pp./filter_TF;
Pvk   = Pv.*scan.w;
ke=fe/scan.w;
Psk = Pvk.*(2*pi*ke).^2;
kmax=fpump./scan.w;
[epsilon,kc]=eps1_mmp(ke,Psk,scan.kvis,kmax);
fc=kc.*scan.w;