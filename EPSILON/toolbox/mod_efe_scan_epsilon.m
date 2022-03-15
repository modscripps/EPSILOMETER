function [Ps_volt_f,Ps_shear_k,Ps_shear_co_k,epsilon,epsilon_co,f,k,fc,kc,Ppan,Ppan_co,fom]=mod_efe_scan_epsilon(scan,shear_channel,Meta_Data)
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
Meta_Data.PROCESS.dof_coh=10;
dof_coh=Meta_Data.PROCESS.dof_coh;
nb_vibration_channel=3;% hard coded nb of acceleration channel used to mv correction


%ALB little warning about dof and statistical significance
% following Lueck2021c: bias in coherent noise correction.
dof=2*length(scan.(shear_channel))/nfft-1;

% If there is data in this shear channel, compute epsilon
if isfinite(scan.(shear_channel))
    
    % -------------------------------------------------------------------------
    % Get constants and transfer functions
    G = 9.81;
    if isfield(Meta_Data,'AFE')
        Sv = Meta_Data.AFE.(shear_channel(1:2)).cal;
    else
        try
            Sv = Meta_Data.epsi.(shear_channel(1:2)).Sv;
        catch
            Sv = Meta_Data.epsi.(shear_channel(1:2)).cal;
        end
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
    Ps_velocity_f = ((2*G/(Sv*w))^2).*Ps_volt_f./filter_TF;
    
    % Convert the frequency velocity spectrum to shear wavenumber spectrum
    %Ps_shear_k = ((2*pi*k).^2).*(Ps_velocity_f./w);
    Ps_shear_k = ((2*pi*k).^2).*(Ps_velocity_f.*w); %NC 9/2/21 - frequency spectrum should be MULTIPLIED by w, not divided
    
    % Compute epsilon using eps1_mmp.m with kmax
    try
        [epsilon,kc(1)]=eps1_mmp(k,Ps_shear_k,scan.kvis,kmax);
        fc(1)=kc(1).*w;
    catch
        epsilon=nan;
        kc(1)=nan;
        fc(1)=nan;
    end
    
    % Remove the coherent part of the frequency spectrum
    
    Ps_velocity_co_f = Ps_velocity_f.*(1-Csa);
    
    if isfield(Meta_Data.PROCESS, "multivariate")
        if (Meta_Data.PROCESS.multivariate)
            % ---------------------------------------------------------------------
            % Now, do the same calculations with the multivariate coherence correction
            % in this mv context 
            % Ps_volt_f_co=Ps_volt_f-Csa
            % Pv_co=(Pv-Coh)./bias;
            Ps_velocity_co_f = Ps_velocity_f.*(1-Csa);
%             bias=1-1.02*nb_vibration_channel./dof_coh;
%             df=Fs./nfft;
%             Ps_volt_co_f=(Ps_volt_f(:)-Csa(:).*df)./bias;
%             Ps_volt_co_f(Ps_volt_co_f<0)=nan;
%             Ps_volt_co_f=fillmissing(Ps_volt_co_f,'linear');
%             close all;
%             Ps_velocity_co_f = ((2*G/(Sv*w))^2).*Ps_volt_co_f./filter_TF;

        end
    end

    
    % Convert the frequency velocity spectrum to shear wavenumber spectrum
    %Ps_shear_co_k = ((2*pi*k).^2).*(Ps_shear_co_f./w);
    
    Ps_shear_co_k = ((2*pi*k).^2).*(Ps_velocity_co_f.*w); %NC 9/2/21 - frequency spectrum should be MULTIPLIED by w, not divided
    
    % Compute epsilon using eps1_mmp.m with kmax
    try
        [epsilon_co,kc(2)]=eps1_mmp(k,Ps_shear_co_k,scan.kvis,kmax);
        fc(2)=kc(2).*w;
    catch
        epsilon_co=nan;
        kc(2)=nan;
        fc(2)=nan;
    end
    
            % Get Panchev spectrum
        if ~isempty(epsilon) && ~isnan(epsilon)
            
            sig_lnS=5/4*dof^(-7/9);
            [kpan,Ppan] = panchev(epsilon,scan.kvis);
            [~,Ppan_co] = panchev(epsilon_co,scan.kvis);
            Ppan_co=interp1(kpan(~isnan(Ppan_co)),Ppan_co(~isnan(Ppan_co)),k);
            Ppan=interp1(kpan(~isnan(Ppan)),Ppan(~isnan(Ppan)),k);

            fom=log(Ps_shear_co_k./Ppan_co);
            fom=nanvar(fom(k>2 & k<kc(2)))./sig_lnS;
%             close all
%             loglog(k,Ps_shear_co_k,'b',k,Ppan,'r')
%             pause
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