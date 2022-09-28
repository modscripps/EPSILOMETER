function scan=mod_epsilometer_make_scan_v2(Profile,scan,Meta_Data)

tscan=Meta_Data.PROCESS.tscan;
Fs_epsi=Meta_Data.PROCESS.Fs_epsi;
Fs_ctd=Meta_Data.PROCESS.Fs_ctd;
channels=Meta_Data.PROCESS.channels;
G=9.91;
twoG=2*G;

% L_epsi=length(Profile.t1);% length of profile
% L_ctd=length(Profile.T);% length of profile
% numbuer of samples for a scan. I make sure it is always even
N_epsi=tscan.*Fs_epsi-mod(tscan*Fs_epsi,2);
N_ctd=tscan.*Fs_ctd-mod(tscan*Fs_ctd,2);

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

[~,indP] = sort(abs(Profile.P-scan.Pr));
indP=indP(1);

ind_Pr_epsi = find(Profile.epsitime<Profile.ctdtime(indP),1,'last');
ind_scan = ind_Pr_epsi-N_epsi/2:ind_Pr_epsi+N_epsi/2; % ind_scan is even

if length(ind_scan)<975
fprintf('length indscan : %i\n',length(ind_scan))
end

% compute mean values of w,T,S of each scans
% scan.w    = average_scan(pr_w,ind_scan);
% scan.t    = average_scan(Profile.T,ind_scan);
% scan.s    = average_scan(Profile.S,ind_scan);
% scan.pr   = average_scan(Profile.P,ind_scan); % I know this redondant with Pr axis
% scan.dnum = average_scan(Profile.time,ind_scan);
% scan.kvis=nu(scan.s,scan.t,scan.pr);
% scan.ktemp=kt(scan.s,scan.t,scan.pr);
% scan.dP=Profile.P(ind_scan([1 end]));

% scan.nfft=Meta_Data.PROCESS.nfft;
% scan.nfftc=Meta_Data.PROCESS.nfftc;
% scan.N=N_epsi;

%compute spectra for acceleration channels.
for c=[inda1 inda2 inda3]
    wh_channels=channels{c};
    scan.(wh_channels)=Profile.(wh_channels)(ind_scan)*G; % time series in m.s^{-2}
end
for c=1:length(channels)
    wh_channel=channels{c};
    switch wh_channel
        case 't1'
            scan.(wh_channel)=Profile.(wh_channel)(ind_scan).*dTdV(1); % time series in Celsius
        case 't2'
            scan.(wh_channel)=Profile.(wh_channel)(ind_scan).*dTdV(2); % time series in Celsius
        case 's1'
            scan.(wh_channel)=Profile.(wh_channel)(ind_scan).*twoG./(Sv1.*scan.w); % time series in m.s^{-1}
        case 's2'
            scan.(wh_channel)=Profile.(wh_channel)(ind_scan).*twoG ./(Sv2.*scan.w); % time series in m.s^{-1}
    end
end
