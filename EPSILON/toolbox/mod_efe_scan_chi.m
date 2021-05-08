function [Pt_volt_f,Pt_Tg_k,chi,f,k,fc,kc,flag_tg_fc]=mod_efe_scan_chi(scan,fpo7_channel,Meta_Data,h_freq,FPO7noise)


% -------------------------------------------------------------------------
% Get constants and transfer functions

w = scan.w;
nfft=Meta_Data.PROCESS.nfft;
try
Fs=Meta_Data.PROCESS.Fs_epsi;
catch
   Fs=Meta_Data.AFE.FS; 
end

switch fpo7_channel
    case 't1_volt'
        if isfield(Meta_Data,'AFE')
            dTdV = Meta_Data.AFE.t1.cal;
        else
        dTdV=Meta_Data.epsi.t1.dTdV;
        end
    case 't2_volt'
                if isfield(Meta_Data,'AFE')
            dTdV = Meta_Data.AFE.t2.cal;
        else
        dTdV=Meta_Data.epsi.t2.dTdV;
                end
    otherwise
        disp('wrong channel to compute chi, must be t1 or t2')
end

% If FPO7 noise is not specified, get it from Meta_Data
if nargin<5
    h_freq=get_filters_MADRE(Meta_Data,f);
    % get FPO7 channel average noise to compute chi
    switch Meta_Data.MAP.temperature
        case 'Tdiff'
            FPO7noise=load(fullfile(Meta_Data.CALIpath,'FPO7_noise.mat'),'n0','n1','n2','n3');
        otherwise
            FPO7noise=load(fullfile(Meta_Data.CALIpath,'FPO7_notdiffnoise.mat'),'n0','n1','n2','n3');
    end
end

filter_TF=h_freq.FPO7(scan.w);

% ---------------------------------------------------------------------
% Calculate spectra and chi

% Compute the frequency spectrum of timeseries in volts
[Pt_volt_f,f] = pwelch(detrend(scan.(fpo7_channel)),nfft,[],nfft,Fs,'psd');
k = f./w;

% Convert frequency spectrum of volt timeseries to frequency spectrum of
% temperature in C
Pt_T_f = (Pt_volt_f*(dTdV^2)) ./ filter_TF;

% Convert temperature frequency spectrum to temperature gradient wavenumber spectrum
Pt_Tg_k = ((2*pi*k).^2).*Pt_T_f./w;

% Calculate chi
dk = nanmean(diff(k));
fc_index = FPO7_cutoff(f,Pt_volt_f,FPO7noise);
fc = f(find(f<=f(fc_index),1,'last'));
kc = fc/w;
krange = find(k<=k(fc_index));
chi = 6*scan.ktemp*dk.*nansum(Pt_Tg_k(krange));

% high signal flag: The cut off frequency is very high. 
% this could mean that the whole scan is corrupt since the spectrum is way
% above the noise floor.
flag_tg_fc = fc_index<round(.95*length(f));
