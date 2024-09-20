function [Cu1a,Cu2a,sumCu1a,sumCu2a,fe]=mod_efe_scan_coherence(scan,acceleration_channel,Meta_Data)
% Compute the coherence over the whole profile
% over the 1./tsan:Fs frequency frequency axis with nfft samples.

nfft=Meta_Data.PROCESS.nfft;

% NC - adding a check if timeseries is at least nfft long
% ALB t1 is not always here (e.g.,NISKINE2018)
% Fnames=fieldnames(scan);
wh_field='ind_scan';
if length(scan.(wh_field))<=nfft
    Cu1a = [];
    Cu2a = [];
    sumCu1a = [];
    sumCu2a = [];
    fe = [];
elseif length(scan.(wh_field))>nfft
    
    try
        Fs=Meta_Data.PROCESS.Fs_epsi;
    catch
        Fs=Meta_Data.AFE.FS;
    end
    fc1=Meta_Data.PROCESS.fc1;
    fc2=Meta_Data.PROCESS.fc2;
    
    if isfinite(scan.s1_volt)
        [Cu1a,fe] = mscohere(detrend(scan.s1_volt),detrend(scan.(acceleration_channel)),nfft,[],nfft,Fs);
        sumCu1a=sum(Cu1a(fe>fc1 & fe<fc2))*mean(diff(fe),'omitmissing');
    else
        fe = nan(nfft/2 + 1,1);
        Cu1a = nan(nfft/2 + 1,1);
        sumCu1a = nan;
    end
    
    if isfinite(scan.s2_volt)
        [Cu2a,~] = mscohere(detrend(scan.s2_volt),detrend(scan.(acceleration_channel)),nfft,[],nfft,Fs);
        sumCu2a=sum(Cu2a(fe>fc1 & fe<fc2))*mean(diff(fe),'omitmissing');
    else
        Cu2a = nan(nfft/2 + 1,1);
        sumCu2a = nan;
    end
    
end %NC - end check that t1_volt is at least nfft long