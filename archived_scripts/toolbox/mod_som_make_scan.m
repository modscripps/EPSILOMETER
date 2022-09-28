function scan=mod_som_make_scan(CTDProfile,EPSIProfile,Pr,Meta_Data)

average_scan=@(x,y) (nanmean(x(y)));


tscan=Meta_Data.tscan;
df_epsi=Meta_Data.df_epsi;
df_ctd=Meta_Data.df_ctd;
channels=Meta_Data.PROCESS.channels;
G=9.91;

switch Meta_Data.vehicle
    case 'FISH'
        pr_w = compute_fallrate_downcast(CTDProfile);
    case 'WW'
        % TODO: make the P from the WW CTD in the same unit as SEABIRD
        pr_w = compute_speed_upcast(CTDProfile);
        pr_w = -pr_w/1e7;
end

LCTD=length(CTDProfile.P);% length of profile
LEPSI=length(EPSIProfile.t1);% length of profile
% numbuer of samples for a scan. I make sure it is always even
N_epsi=tscan.*df_epsi-mod(tscan*df_epsi,2);
N_ctd=tscan.*df_ctd-mod(tscan*df_ctd,2);

% get index of the acceleration channels (usefull when the number of channels is not 8)
inda1=find(cellfun(@(x) strcmp(x,'a1'),channels));
inda2=find(cellfun(@(x) strcmp(x,'a2'),channels));
inda3=find(cellfun(@(x) strcmp(x,'a3'),channels));

% get shear probe calibration coefficient.
Sv1=Meta_Data.epsi.s1.Sv;
Sv2=Meta_Data.epsi.s2.Sv;
% Sensitivity of FPO7 probe, nominal
dTdV(1)=Meta_Data.epsi.t1.dTdV; % define in mod_epsi_temperature_spectra
dTdV(2)=Meta_Data.epsi.t2.dTdV; % define in mod_epsi_temperature_spectra 


%% make the scan

[~,indP] = sort(abs(CTDProfile.P-Pr));
indP=indP(1);
ind_ctdscan = indP-N_ctd/2:indP+N_ctd/2; % ind_scan is even

ind_Pr_epsi = find(EPSIProfile.epsitime<CTDProfile.ctdtime(indP),1,'last');
ind_epsiscan = ind_Pr_epsi-N_epsi/2:ind_Pr_epsi+N_epsi/2; % ind_scan is even

ind_ctdscan=ind_ctdscan(ind_ctdscan>0 & ind_ctdscan<LCTD);
ind_epsiscan=ind_epsiscan(ind_epsiscan>0 & ind_epsiscan<LEPSI);

scan.w   = average_scan(pr_w,ind_ctdscan(ind_ctdscan>0 & ind_ctdscan<LCTD));

% compute mean values of w,T,S of each scans
scan.w    = average_scan(pr_w,ind_ctdscan);
scan.t    = average_scan(CTDProfile.T,ind_ctdscan);
scan.s    = average_scan(CTDProfile.S,ind_ctdscan);
scan.pr   = average_scan(CTDProfile.P,ind_ctdscan); % I know this redondant with Pr axis
scan.dnum = average_scan(CTDProfile.ctdtime,ind_ctdscan);
scan.kvis=nu(scan.s,scan.t,scan.pr);
scan.ktemp=kt(scan.s,scan.t,scan.pr);
scan.dP=CTDProfile.P(ind_ctdscan([1 end]));

scan.nfft=Meta_Data.PROCESS.nfft;
scan.nfftc=Meta_Data.PROCESS.nfftc;
scan.N=N_epsi;

%compute spectra for acceleration channels.
for c=[inda1 inda2 inda3]
    wh_channels=channels{c};
    scan.(wh_channels)=EPSIProfile.(wh_channels)(ind_epsiscan)*G; % time series in m.s^{-2}
end
for c=1:length(channels)
    wh_channel=channels{c};
    switch wh_channel
        case 't1'
            scan.(wh_channel)=EPSIProfile.(wh_channel)(ind_epsiscan).*dTdV(1); % time series in Celsius
        case 't2'
            scan.(wh_channel)=EPSIProfile.(wh_channel)(ind_epsiscan).*dTdV(2); % time series in Celsius
        case 's1'
            scan.(wh_channel)=EPSIProfile.(wh_channel)(ind_epsiscan).*G./(Sv1.*scan.w); % time series in m.s^{-1}
        case 's2'
            scan.(wh_channel)=EPSIProfile.(wh_channel)(ind_epsiscan).*G ./(Sv2.*scan.w); % time series in m.s^{-1}
    end
end

