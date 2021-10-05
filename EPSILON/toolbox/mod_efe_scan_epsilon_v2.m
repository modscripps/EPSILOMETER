function [Ps_volt_f,Ps_shear_k,Ps_shear_co_k,epsilon,epsilon_co,f,k,fc,kc,Ppan,Ppan_co,fom]=mod_efe_scan_epsilon_v2(scan,shear_channel,Meta_Data)
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

nfft=Meta_Data.PROCESS.nfft;
% create a matrix with the accel channel we want to use to correct the
% shear signal
% vibration=[scan.a3_g scan.a1_g scan.a2_g];
vibration=scan.a3_g;
nb_vib_channel=size(vibration,2);

%ALB little warning about dof and statistical significance
% following Lueck2021c: bias in coherent noise correction.
dof=2*length(scan.(shear_channel))/nfft-1;
if dof<3+nb_vib_channel
    disp(["dof is too small. Please consider reducing the number of " ...
            "accell channels for the vibration correction or increase tscan"]);
end
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
    
    if isfield(Meta_Data,'AFE')
        Fs=Meta_Data.AFE.FS;
    else
        Fs=Meta_Data.PROCESS.Fs_epsi;
    end
    fpump=Meta_Data.PROCESS.ctd_fc;
    kmax=fpump./w;
    
    % Get the filter transfer functions.
    h_freq=Meta_Data.PROCESS.h_freq;
    
    
    % ---------------------------------------------------------------------
    % Calculate spectra and epsilon
    
    % Compute the frequency spectrum of timeseries in volts
%     [Ps_volt_f,f] = pwelch(detrend(scan.(shear_channel)),nfft,[],nfft,Fs,'psd');
    [Ps_volt_f,Ps_volt_co_f,f]= ...
                 multivariate_correction(scan.(shear_channel),vibration,nfft,Fs);
    
    % multivariate correction of shear signal using all three
    % accelerometers channels.
    
    k = f./w;
    filter_TF=(h_freq.shear .* haf_oakey(f,w));
    
    % Convert frequency spectrum of velocity timeseries in volts to the frequency spectrum
    % of shear in s^-1
    Ps_velocity_f     = ((2*G/(Sv*w))^2).*Ps_volt_f./filter_TF;
    Ps_velocity_co_f  = ((2*G/(Sv*w))^2).*Ps_volt_co_f./filter_TF;
    
    % Convert the frequency velocity spectrum to shear wavenumber spectrum
    Ps_shear_k    = ((2*pi*k).^2).*(Ps_velocity_f.*w);    
    Ps_shear_co_k = ((2*pi*k).^2).*(Ps_velocity_co_f.*w); 

    % Compute epsilon using eps1_mmp.m with kmax
    try
        [epsilon,kc(1)]=eps1_mmp(k,Ps_shear_k,scan.kvis,kmax);
        fc(1)=kc(1).*w;
    catch
        epsilon=nan;
        kc(1)=nan;
        fc(1)=nan;
    end
    % Compute epsilon using eps1_mmp.m with kmax
    try
        [epsilon_co,kc(1)]=eps1_mmp(k,Ps_shear_co_k,scan.kvis,kmax);
        fc(1)=kc(1).*w;
    catch
        epsilon_co=nan;
        kc(1)=nan;
        fc(1)=nan;
    end
    
                % Get Panchev spectrum
        if ~isempty(epsilon)
            
            sig_lnS=5/4*dof^(-7/9);
            [kpan,Ppan] = panchev(epsilon,scan.kvis);
            [~,Ppan_co] = panchev(epsilon_co,scan.kvis);
            Ppan_co=interp1(kpan(~isnan(Ppan_co)),Ppan_co(~isnan(Ppan_co)),k);
            Ppan=interp1(kpan(~isnan(Ppan)),Ppan(~isnan(Ppan)),k);

            fom=log(Ps_shear_co_k./Ppan_co);
            fom=nanvar(fom(k<kc(1)))./sig_lnS;
        else
            fom=NaN;
            Ppan_co=NaN;
            Ppan=NaN;
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
    fom=NaN;
    Ppan_co=NaN;
    Ppan=NaN;

    
end