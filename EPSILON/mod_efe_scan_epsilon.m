% function [Ps_volt_f,Ps_velocity_f,Ps_shear_k,Ps_shear_co_k,...
%           epsilon,epsilon_co,method,f,k,fc,kc, ...
%           Ppan,Ppan_co,fom,calib_volt,calib_vel]= ...
%           mod_efe_scan_epsilon(scan,shear_channel,Meta_Data)
function scan=mod_efe_scan_epsilon(scan,shear_channel,Meta_Data)

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
% Meta_Data.PROCESS.dof_coh=10;
% dof_coh=Meta_Data.PROCESS.dof_coh;
% nb_vibration_channel=3;% hard coded nb of acceleration channel used to mv correction


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

    try
        Fs=Meta_Data.AFE.FS;
    catch
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
            Csa    = scan.Cs1a.a3;
        case 's2_volt'
            Csa    = scan.Cs2a.a3;
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


    % Another qc flag: ratio between shear./a3 at the CTD peak.
    % Get idx with fpump - 5Hz < f < fpump + 5Hz
    % find the max FOCO for shear and a3 within this range and compute the
    % ratio. Keep the ration


    idx=find(scan.f>fpump-5 & scan.f<fpump+5);
    calib_volt=max(Ps_volt_f(idx))./max(scan.Pa_g_f.a3(idx));
    calib_vel=max(Ps_velocity_f(idx))./max(scan.Pa_g_f.a3(idx));

%     if abs(calib_vel-0.4)>0.4
%         Ps_velocity_f=Ps_velocity_f .* 0.5/calib_vel;
%     end

    % Convert the frequency velocity spectrum to shear wavenumber spectrum
    Ps_shear_k = ((2*pi*k).^2).*(Ps_velocity_f.*w);

    % Compute epsilon using eps1_mmp.m with kmax
    try
        [epsilon,kmin(1),kc(1),method(1)]=eps1_mmp(k,Ps_shear_k,scan.kvis,kmax);
        fc(1)=kc(1).*w;
        if isinf(epsilon)
            epsilon=nan;
        end

    catch
        epsilon=nan;
        kc(1)=nan;
        fc(1)=nan;
        kmin(1)=nan;
        method(1)=nan;
    end

    % Remove the coherent part of the frequency spectrum
    if isfield(scan,'Cs1a')
        Ps_velocity_co_f = Ps_velocity_f.*(1-Csa);
    end


    % Convert the frequency velocity spectrum to shear wavenumber spectrum
    %Ps_shear_co_k = ((2*pi*k).^2).*(Ps_shear_co_f./w);

    if isfield(scan,'Cs1a')
        Ps_shear_co_k = ((2*pi*k).^2).*(Ps_velocity_co_f.*w);


        % Compute epsilon using eps1_mmp.m with kmax
        try
            [epsilon_co,kmin(2),kc(2),method(2)]=eps1_mmp(k,Ps_shear_co_k,scan.kvis,kmax);
            if isinf(epsilon_co)
                epsilon_co=nan;
            end
            fc(2)=kc(2).*w;
        catch
            epsilon_co=nan;
            kc(2)=nan;
            fc(2)=nan;
            kmin(2)=nan;
            method(2)=nan;
        end

    else
        Ps_shear_co_k = nan;
        epsilon_co = nan;
        kc(2) = nan;
        fc(2) = nan;
        kmin(2)=nan;
        method(2)=nan;
    end

    % Get Panchev spectrum
    if ~isempty(epsilon) && ~isnan(epsilon)

        % [kpan,Ppan] = panchev(epsilon,scan.kvis);
        % Ppan=interp1(kpan(~isnan(Ppan)),Ppan(~isnan(Ppan)),k);

        Pxx=Ps_shear_co_k(k>=kmin(2) & k<=kc(2));
        kin=k(k>=kmin(2) & k<=kc(2));

        sig_lnS=sqrt(5/4*(dof-1)^(-7/9));
        Ns=2*length(k);
        Tm=0.8+sqrt(1.56/Ns);
        [kpan,Ppan_co] = panchev(epsilon_co,scan.kvis);
        interp_Ppan=interp1(kpan(~isnan(Ppan_co)),Ppan_co(~isnan(Ppan_co)),kin);
        % Ppan_co=interp1(kpan(~isnan(Ppan_co)),Ppan_co(~isnan(Ppan_co)),k);
        
        fom=log10(Pxx(:)./interp_Ppan(:));
        mad_spec=mad(fom);
        fom=mad_spec./sig_lnS./Tm;

    else
        fom=NaN;
        % Ppan_co=NaN;
        % Ppan=NaN;
    end

else

    Ps_volt_f = nan(nfft/2 + 1,1);
    Ps_velocity_f = nan(nfft/2 + 1,1);
    Ps_shear_k = nan(nfft/2 + 1,1);
    Ps_shear_co_k = nan(nfft/2 + 1,1);

    epsilon = nan;
    epsilon_co = nan;
    method = nan;
    fc = nan(nfft/2 + 1,1);
    kc = nan;
    kmin = nan;
    fom=NaN;
    % Ppan_co=NaN;
    % Ppan=NaN;
    calib_volt=nan;
    calib_vel=nan;

end

% % Put new variables in the structure
% varList = {'Ps_volt_f','Ps_velocity_f','Ps_shear_k','Ps_shear_co_k',...
%     'epsilon','epsilon_co','fc','kc','Ppan','Ppan_co','fom','calib_volt','calib_vel','method'};
% for iVar=1:numel(varList)
%     scan.(varList{iVar}).(currChannel(1:2)) = eval(varList{iVar});
% end

scan.Ps_volt_f.(shear_channel(1:2))     = Ps_volt_f;
scan.Ps_velocity_f.(shear_channel(1:2)) = Ps_velocity_f;
scan.Ps_shear_k.(shear_channel(1:2))    = Ps_shear_k;
scan.Ps_shear_co_k.(shear_channel(1:2)) = Ps_shear_co_k;
scan.epsilon.(shear_channel(1:2))       = epsilon;
scan.epsilon_co.(shear_channel(1:2))    = epsilon_co;
scan.method.(shear_channel(1:2))        = method;
scan.fc.(shear_channel(1:2))            = fc;
scan.kc.(shear_channel(1:2))            = kc;
scan.kmin.(shear_channel(1:2))          = kmin;
scan.fom.(shear_channel(1:2))           = fom;
scan.calib_volt.(shear_channel(1:2))    = calib_volt;
scan.calib_vel.(shear_channel(1:2))     = calib_vel;

