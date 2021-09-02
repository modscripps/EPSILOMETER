function [Ps_volt_f,Ps_shear_k,Ps_shear_co_k,epsilon,epsilon_co,f,k,fc,kc]=mod_efe_scan_epsilon(scan,shear_channel,acceleration_channel,Meta_Data)
% get epsilon and the cutoff frequency
%
% OUTPUTS:
%   P   = shear frequency power spectrum
%   Pv  = non-coherent shear frequency power spectrum (full profile coherence with a3 channel has been removed)
%   Pvk = non-coherent shear wavenumber power spectrum (not saved)
%   Psk = non-coherent shear 2pi*wavenumber power spectrum
%   Csa = full profile coherence between shear channel and a3, computed
%         earlier with mod_efe_scan_coherence
%   epsilon    = epsilon calculated from Psk
%   epsilon_co = epsilon calculated from Psk
%   fc  = cutoff frequency
%   f  = frequency array

% If there is data in this shear channel, compute epsilon
if isfinite(scan.(shear_channel))
    
    % -------------------------------------------------------------------------
    % Get constants and transfer functions
    G = 9.81;
    if isfield(Meta_Data,'AFE')
        Sv = Meta_Data.AFE.(shear_channel(1:2)).cal;
    else
        Sv = Meta_Data.epsi.(shear_channel(1:2)).Sv;
    end
    % NC - Changed w to absolute value so this also works for upcasts
    w = abs(scan.w);
    nfft=Meta_Data.PROCESS.nfft;
    if isfield(Meta_Data,'AFE')
        Fs=Meta_Data.AFE.FS;
    else
        Fs=Meta_Data.PROCESS.Fs_epsi;
    end
    fpump=Meta_Data.PROCESS.ctd_fc;
    kmax=fpump./w;
    
    % Get the filter transfer functions.
    h_freq=Meta_Data.PROCESS.h_freq;
    
    % Get the coherence between the shear channel and acceleration channel (should
    % be a3)
    switch shear_channel
        case 's1_volt'
            Csa = scan.Cs1a3_full;
        case 's2_volt'
            Csa = scan.Cs2a3_full;
    end
    
    % ---------------------------------------------------------------------
    % Calculate spectra and epsilon
    
    % Compute the frequency spectrum of timeseries in volts
    [Ps_volt_f,f] = pwelch(detrend(scan.(shear_channel)),nfft,[],nfft,Fs,'psd');
    k = f./w;
    filter_TF=(h_freq.shear .* haf_oakey(f,w));
    
    % Convert frequency spectrum of velocity timeseries in volts to the frequency spectrum
    % of shear in s^-1
    Ps_shear_f = ((2*G/(Sv*w))^2).*Ps_volt_f./filter_TF;
    
    % Convert the shear frequency spectrum to shear wavenumber spectrum
    %Ps_shear_k = ((2*pi*k).^2).*(Ps_shear_f./w);
    Ps_shear_k = ((2*pi*k).^2).*(Ps_shear_f.*w); %NC 9/2/21 - frequency spectrum should be MULTIPLIED by w, not divided
    
    % Compute epsilon using eps1_mmp.m with kmax
    try
        [epsilon,kc(1)]=eps1_mmp(k,Ps_shear_k,scan.kvis,kmax);
        fc(1)=kc(1).*w;
    catch
        epsilon=nan;
        kc(1)=nan;
        fc(1)=nan;
    end
    
    
    % ---------------------------------------------------------------------
    % Now, do the same calculations with the coherence correction
    
    % Remove the coherent part of the frequency spectrum
    Ps_shear_co_f = Ps_shear_f.*(1-Csa);
    
    % Convert the frequency velocity spectrum to shear wavenumber spectrum
    %Ps_shear_co_k = ((2*pi*k).^2).*(Ps_shear_co_f./w);
    Ps_shear_co_k = ((2*pi*k).^2).*(Ps_shear_co_f.*w); %NC 9/2/21 - frequency spectrum should be MULTIPLIED by w, not divided
    
    % Compute epsilon using eps1_mmp.m with kmax
    try
        [epsilon_co,kc(2)]=eps1_mmp(k,Ps_shear_co_k,scan.kvis,kmax);
        fc(1)=kc(2).*w;
    catch
        epsilon_co=nan;
        kc(1)=nan;
        fc(1)=nan;
    end
    
    
else
    
    Ps_volt_f = nan(nfft/2 + 1,1);
    Ps_shear_k = nan(nfft/2 + 1,1);
    Ps_shear_co_k = nan(nfft/2 + 1,1);
    epsilon = nan;
    epsilon_co = nan;
    f = nan(nfft/2 + 1,1);
    k = nan;
    fc = nan(nfft/2 + 1,1);
    kc = nan;
    
end