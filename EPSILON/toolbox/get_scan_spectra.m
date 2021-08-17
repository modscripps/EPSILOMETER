function [scan] = get_scan_spectra(Profile,id_scan)

% function [scan] = get_scan_spectra(Profile,id_scan)
%
% INPUTS:
%   Profile = Profile### structure
%   id_scan = Index of Profile.pr on which to compute spectra
%
% OUTPUTS:
%   scan = Structure similar to Profile### but only for current index and
%          now including spectra
%
% Nicole Couto | June 2020
% -------------------------------------------------------------------------
%% Gather some data
% -------------------------------------------------------------------------
Meta_Data   = Profile.Meta_Data;

Pr          = Profile.pr(id_scan);

fpump       = Meta_Data.PROCESS.ctd_fc;
tscan       = Meta_Data.PROCESS.tscan;
try
    Fs_epsi = Meta_Data.PROCESS.Fs_epsi;
catch
    Fs_epsi = Meta_Data.AFE.FS;
end
N_epsi      = tscan.*Fs_epsi-mod(tscan*Fs_epsi,2);
Fs_ctd      = Meta_Data.PROCESS.Fs_ctd;
N_ctd       = tscan.*Fs_ctd-mod(tscan*Fs_ctd,2);
h_freq      = Meta_Data.PROCESS.h_freq;

% Find indices of ctd and epsi data that make up this scan
[~,indP]    = sort(abs(Profile.ctd.P-Pr));
indP        = indP(1);
ind_ctdscan = indP-N_ctd/2:indP+N_ctd/2; % ind_scan is even
ind_Pr_epsi = find(Profile.epsi.time_s<Profile.ctd.time_s(indP),1,'last');
ind_scan    = ind_Pr_epsi-N_epsi/2:ind_Pr_epsi+N_epsi/2; % ind_scan is even


% get FPO7 channel average noise to compute chi
try 
    tempChoice = Meta_Data.MAP.temperature;
catch
    tempChoice = Meta_Data.AFE.temp_circuit;
end

if Meta_Data.PROCESS.adjustTemp
    
    FPO7noise = Meta_Data.MAP.Tnoise1;
    
else
    
    switch tempChoice
        case 'Tdiff'
            Meta_Data.PROCESS.FPO7noise=load(fullfile(Meta_Data.CALIpath,'FPO7_noise.mat'),'n0','n1','n2','n3');
        otherwise
            Meta_Data.PROCESS.FPO7noise=load(fullfile(Meta_Data.CALIpath,'FPO7_notdiffnoise.mat'),'n0','n1','n2','n3');
    end
    FPO7noise   = Meta_Data.PROCESS.FPO7noise;
    
end

channels    = Meta_Data.PROCESS.timeseries;
LCTD        = length(Profile.ctd.P);% length of profile
% NC 17 July 2021 - scan.w should really be in m/s. I added steps in
% mod_som_read_epsi_files_v3.m to convert pressure to depth (z) and
% calculate dzdt.  Now I'm using dzdt instead of dPdt to deefine scan.w.
%scan.w      = nanmean(Profile.ctd.dPdt(ind_ctdscan(ind_ctdscan>0 & ind_ctdscan<LCTD)));
if isfield(Profile,'dzdt')
    scan.w      = nanmean(Profile.ctd.dzdt(ind_ctdscan(ind_ctdscan>0 & ind_ctdscan<LCTD)));
else
    scan.w      = nanmean(Profile.ctd.dPdt(ind_ctdscan(ind_ctdscan>0 & ind_ctdscan<LCTD)));
end

% check if the scan is not too shallow or too close to the end of the
% profile.
%
% NC - removing the minimum speed check. Sometimes we want up profiles!
% Sometimes we want both. If we add this back later, I think it should be
% coded into Meta_Data, not hard-coded here. I will leave a criteria that
% scan.w needs to be not nan and not inf.
%
% limit_speed = 0.2;
% if ind_ctdscan(1)>0 && ind_ctdscan(end)<=length(Profile.ctd.time_s) && scan.w>limit_speed ...
%         && ind_scan(1)>0 && ind_scan(end)<=length(Profile.epsi.time_s)
if ind_ctdscan(1)>0 && ind_ctdscan(end)<=length(Profile.ctd.time_s) ...
        && ind_scan(1)>0 && ind_scan(end)<=length(Profile.epsi.time_s) ...
        && ~isinf(scan.w) && ~isnan(scan.w);
    
    % Put new variables in the structure
    varList = {'Pr','tscan','Fs_epsi','N_epsi',...
        'Fs_ctd','N_ctd','h_freq',...
        'indP','ind_ctdscan','ind_scan','FPO7noise'};
    for iVar=1:numel(varList)
        scan.(varList{iVar}) = eval(varList{iVar});
    end
    
    scan.pr     = nanmean(Profile.ctd.P(ind_ctdscan));
    scan.t      = nanmean(Profile.ctd.T(ind_ctdscan));
    scan.s      = nanmean(Profile.ctd.S(ind_ctdscan));
    if isfield(Profile.ctd,'dnum')
        scan.dnum   = nanmean(Profile.ctd.dnum(ind_ctdscan));
    else
        scan.time_s = nanmean(Profile.ctd.time_s(ind_ctdscan));
    end
    % NC 17 July 2021 - the functions nu and kt want pressure in MPa, not
    % db. Convert to MPa before calculating kinematic viscosity and thermal
    % diffusivity.
%     scan.kvis   = nu(scan.s,scan.t,scan.pr);
%     scan.ktemp  = kt(scan.s,scan.t,scan.pr);
    scan.kvis   = nu(scan.s,scan.t,db2MPa(scan.pr));
    scan.ktemp  = kt(scan.s,scan.t,db2MPa(scan.pr));
    scan.kmax   = fpump./scan.w; 
    
    % Add timeseries for each channel 
    for c=1:length(channels)
        currChannel=channels{c};
        scan.(currChannel)=Profile.epsi.(currChannel)(ind_scan); % time series in m.s^{-2}
    end
    
    % Add full profile coherences
    scan.Cs1a3_full = Profile.Cs1a3_full;
    scan.Cs2a3_full = Profile.Cs2a3_full; 
    
    %% Compute spectra for acceleration channels
    % ---------------------------------------------------------------------
    
    idxA = contains(Meta_Data.PROCESS.timeseries,'a');
    chanList = Meta_Data.PROCESS.timeseries(idxA);
    for iChan=1:numel(chanList)
        %c = find(cellfun(@(x) strcmp(x,chanList{iChan}),channels));
        %currChannel = channels{c};
        currChannel = chanList{iChan};
        
        % Get the spectrum for the current acceleration channel
        [Pa_g_f,sumPa_g,f]= mod_efe_scan_acceleration(scan,currChannel,Meta_Data);
        
        % Get the coherence spectra between each shear probe and the current
        % acceleration channel
        [Cs1a,Cs2a,sumCs1a,sumCs2a,~]= ...
            mod_efe_scan_coherence(scan,currChannel,Meta_Data);
        
        % Put new variables in the structure
        varList = {'Pa_g_f','sumPa_g','f','Cs1a','Cs2a','sumCs1a','sumCs2a'};
        for iVar=1:numel(varList)
            scan.(varList{iVar}).(currChannel(1:2)) = eval(varList{iVar});
        end
    end
    
    scan.accelNoise=45e-6^2+0*scan.f.a1;
    
    %% Compute shear spectra and epsilon
    % ---------------------------------------------------------------------
    
    idxS = contains(Meta_Data.PROCESS.timeseries,'s');
    chanList = Meta_Data.PROCESS.timeseries(idxS);
    for iChan=1:numel(chanList)
        currChannel = chanList{iChan};
        
        [Ps_volt_f,Ps_shear_k,Ps_shear_co_k,epsilon,epsilon_co,f,k,fc,kc] = ...
            mod_efe_scan_epsilon(scan,currChannel,'a3',Meta_Data);
        
        % Get Panchev spectrum
        if ~isempty(epsilon)
            [kpan,Ppan] = panchev(epsilon,scan.kvis);
            [~,Ppan_co] = panchev(epsilon_co,scan.kvis);
        end
        
        % Put new variables in the structure
        varList = {'Ps_volt_f','Ps_shear_k','Ps_shear_co_k',...
            'epsilon','epsilon_co','f','k','fc','kc','Ppan','Ppan_co','kpan'};
        for iVar=1:numel(varList)
            scan.(varList{iVar}).(currChannel(1:2)) = eval(varList{iVar});
        end
        
    end
    
    
    
    %% Compute temperature spectra and chi
    % ---------------------------------------------------------------------
    
    idxT = contains(Meta_Data.PROCESS.channels,'t');
    chanList = Meta_Data.PROCESS.timeseries(idxT);
    for iChan=1:numel(chanList)
        currChannel = chanList{iChan};
        chanFieldName = sprintf('%sspectra',currChannel);
        
        [Pt_volt_f,Pt_Tg_k,chi,f,k,fc,kc,flag_tg_fc] = ...
            mod_efe_scan_chi(scan,currChannel,Meta_Data,h_freq,FPO7noise);
        
        % Put new variables in the structure
        varList = {'Pt_volt_f','Pt_Tg_k','chi','f','k','fc','kc','flag_tg_fc'};
        for iVar=1:numel(varList)
            scan.(varList{iVar}).(currChannel(1:2)) = eval(varList{iVar});
        end
    end
    
    %%
    scan.f = scan.f.s1;
    scan.k = scan.k.s1;
    
    %% Add variable descriptions and units
    % ---------------------------------------------------------------------
    scan.varInfo.Pr = {'Epsi pressure','db'};
    scan.varInfo.pr = {'CTD pressure','db'};  
    if isfield(Profile,'dzdt')
        scan.varInfo.z = {'CTD depth','m'};
        scan.varInfo.w = {'fall speed','m s^{-1}'};
    else
        scan.varInfo.w = {'fall speed','db s^{-1}'};
    end
    scan.varInfo.t = {'temperature','C'}; 
    scan.varInfo.s = {'salinity','psu'};   
    scan.varInfo.dnum = {'datenum','Matlab datenum'};   
    scan.varInfo.kvis = {'kinematic viscosity',''};  
    scan.varInfo.ktemp = {'',''}; 
    scan.varInfo.kmax = {'',''};      
    scan.varInfo.tscan = {'length of scan window','s'};
    scan.varInfo.Fs_epsi = {'',''};
    scan.varInfo.N_epsi = {'',''};
    scan.varInfo.Fs_ctd = {'',''};
    scan.varInfo.N_ctd = {'',''};
    scan.varInfo.h_freq = {'',''};
    scan.varInfo.indP = {'',''};   
    scan.varInfo.ind_ctdscan = {'indices of Profile.ctd.time_s in this scan'; ''};
    scan.varInfo.ind_scan  = {'indices of Profile.epsi.time_s in this scan'; ''};
    scan.varInfo.FP07noise = {'',''};
    scan.varInfo.a1_g = {'acceleration sensor 1 timeseries in this scan','[g]'};
    scan.varInfo.a2_g = {'acceleration sensor 2 timeseries in this scan','[g]'};
    scan.varInfo.a3_g = {'acceleration sensor 3 timeseries in this scan','[g]'};
    scan.varInfo.s1_volt = {'shear sensor 1 timeseries in this scan','Volts'};
    scan.varInfo.s2_volt = {'shear sensor 2 timeseries in this scan','Volts'};
    scan.varInfo.t1_volt = {'temperature sensor 1 timeseries in this scan','Volts'};
    scan.varInfo.t2_volt = {'temperature sensor 2 timeseries in this scan','Volts'};
    scan.varInfo.c_count = {'additional sensor timseries in the scan' 'count'};
    scan.varInfo.f = {'frequency','Hz'};
    scan.varInfo.k = {'wavenumber','cycles m^-^1'};
    scan.varInfo.fc  = {'cutoff frequency, 1=uncorrected, 2=coherence-corrected', 'Hz'};
    scan.varInfo.kc = {'cutoff wavenumber, 1=uncorrected, 2=coherence-corrected', 'cpm'};  
    scan.varInfo.Cs1a3_full = {'coherence betwen s1 and a3 channels between Meta_Data.PROCESS.Prmin and Meta_Data.PROCESS.Prmax',''};
    scan.varInfo.Cs2a3_full = {'coherence betwen s2 and a3 channels between Meta_Data.PROCESS.Prmin and Meta_Data.PROCESS.Prmax',''};
    scan.varInfo.Pa_g_f = {'accleration frequency power spectrum', '[g]^2 Hz^{-1}'};
    scan.varInfo.sumPa_g = {'integrated Pa_f between Meta_Data.PROCESS.fc1 and Meta_Data.PROCESS.fc2','[g]^2 Hz^{-1}'};
    scan.varInfo.Cs1a = {'coherence between s1 and acceleration channels','unitless'};
    scan.varInfo.Cs2a = {'coherence between s2 and acceleration channels','unitless'};
    scan.varInfo.sumCs1a = {'integrated Cu1a between Meta_Data.PROCESS.fc1 and Meta_Data.PROCESS.fc2','Volts^2'};
    scan.varInfo.sumCs2a = {'integrated Cu2a between Meta_Data.PROCESS.fc1 and Meta_Data.PROCESS.fc2','Volts^2'};
    scan.varInfo.accelNoise = {'',''};
    scan.varInfo.Ps_volt_f = {'shear frequency power spectrum', 'Volts^2 Hz^{-1}'};
    scan.varInfo.Ps_shear_k = {'shear wavenumber power spectrum', 's{-1} cpm^{-1}'};
    scan.varInfo.Ps_shear_co_k  = {'coherence-corrected shear frequency power spectrum (full profile coherence with a3 channel has been removed)', ''};
    scan.varInfo.epsilon = {'turbulent kinetic energy dissipation rate calculated from Ps_shear_k', ''};
    scan.varInfo.epsilon_co = {'turbulent kinetic energy dissipation rate calculated from Ps_shear_co_k', ''};
    scan.varInfo.Ppan = {'frequency Panchev curve',''};
    scan.varInfo.Ppan_co = {'frequency Panchev curve',''};
    scan.varInfo.kpan = {'wavenumber array for Panchev curve',''};
    scan.varInfo.Pt_volt_f = {'temperature frequency power spectrum','Volts^2 Hz{-1}'};
    scan.varInfo.Pt_Tg_k = {'temperature gradient wavenumber power spectrum', 'C^2 s{-1} cpm^{-1}'};
    scan.varInfo.chi = {'temperature gradient dissipation rate',''};
    scan.varInfo.flag_tg_fc = {'temperature gradient cut off frequency is very high','0/1'};
    
    
end %endif scan is not too shallow or too close to the end of the profile

