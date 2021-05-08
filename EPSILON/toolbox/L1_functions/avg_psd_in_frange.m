function [PSDavg] = avg_psd_in_frange(psd,f,f1,f2)
% [PSDavg] = avg_psd_in_frange(psd,f,f1,f2);
%
% psd = power spectral density [1xn]
% f   = frequency array        [1xn]
% f1  = lower frequency limit
% f2  = upper frequency limit 
%
% Find average values of power spectra or coherence spectra
% between frequencies f1 and f2

idxF = f>=f1 & f<f2;
PSDavg = 10^nanmean(log10(psd(idxF)));
