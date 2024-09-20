function [Pa,sumP,fe]=mod_efe_scan_acceleration(scan,acceleration_channel,Meta_Data)
% get epsilon and the cutting frequency
%
% OUTPUTS
%   Pa      = acceleration frequency power spectrum
%   sumP    = integrated acceleratation frequency power spectrum between
%             Meta_Data.PROCESS.fc1 and Meta_Data.PROCESS.fc2
%   fe      = frequency array

nfft=Meta_Data.PROCESS.nfft;

if isfield(Meta_Data,'AFE') && isclassfield(Meta_Data.AFE,'FS')
    Fs=Meta_Data.AFE.FS;
elseif isclassfield(Meta_Data.PROCESS,'Fs_epsi')
    Fs=Meta_Data.PROCESS.Fs_epsi;
end
fc1=Meta_Data.PROCESS.fc1;
fc2=Meta_Data.PROCESS.fc2;
h_freq=Meta_Data.PROCESS.h_freq;

[P,fe] = pwelch(detrend(scan.(acceleration_channel)),nfft,[],nfft,Fs,'psd');


filter_TF=(h_freq.electAccel);
Pa   = P./filter_TF;
sumP=sum(Pa(fe>fc1 & fe<fc2))*mean(diff(fe),'omitmissing');
