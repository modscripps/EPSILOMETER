function Profile=calc_turbulence(CTDProfile,EpsiProfile,dz,Meta_Data)

%  Profile structure for Micro Structure. Inside Profile you ll find
%  temperature spectra in degC Hz^-1
%  Horizontal  velocity spectra in m^2/s^-2 Hz^-1
%  Acceleration/speed spectra in s^-1 Hz^-1 
%
%  Created by Arnaud Le Boyer on 7/28/18.

%% get channels
channels=Meta_Data.PROCESS.channels;
nb_channels=length(channels);

nfft=Meta_Data.PROCESS.nfft;
nfftc=Meta_Data.PROCESS.nfftc;
df_epsi=Meta_Data.df_epsi;
df_ctd=Meta_Data.df_ctd;


fpump=Meta_Data.fpump;
tscan=Meta_Data.tscan;

%% Gravity  ... of the situation :)
G       = 9.81;
twoG    =2*G;
%% limit speed 20 cm s^{-1}
limit_speed=.2;
%% define the fall rate of the Profile. 
%  Add Profile.w with w the vertical vel. 
%  We are using the pressure from other sensors (CTD);

%% TODO check the dimension of the raw time series to remove this line
Profile=structfun(@(x) x(:),CTDProfile,'un',0);

switch Meta_Data.vehicle
    case 'FISH'
        Profile.dPdt = compute_fallrate_downcast(CTDProfile);
    case 'WW'
        % TODO: make the P from the WW CTD in the same unit as SEABIRD
        Profile.dPdt = compute_speed_upcast(CTDProfile);
        Profile.dPdt=-Profile.dPdt/1e7;
end
%% define a Pressure axis to an which I will compute epsilon and chi. 
%  The spectra will be nfft long centered around P(z) +/- tscan/2. 
%  

Pr=ceil(min(Profile.P)):dz:floor(max(Profile.P));
nbscan=length(Pr);

LCTD=length(CTDProfile.P);% length of profile
LEPSI=length(EpsiProfile.t1);% length of profile
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

% get FPO7 channel average noise to compute chi
switch Meta_Data.MAP.temperature
    case 'Tdiff'
        FPO7noise=load(fullfile(Meta_Data.CALIpath,'FPO7_noise.mat'),'n0','n1','n2','n3');
    otherwise
        FPO7noise=load(fullfile(Meta_Data.CALIpath,'FPO7_notdiffnoise.mat'),'n0','n1','n2','n3');
end



% start creating the Profile structure
% get nb of scans in the profile
Profile.pr        = Pr;
Profile.nbscan    = nbscan;
Profile.nfft      = nfft;
Profile.tscan     = tscan;
Profile.fpump     = fpump; % arbitrary cut off frequency usually extract from coherence spectra shear/accel 
Profile.nbchannel = nb_channels;

%initialize process flags
Profile.process_flag=Pr*0+1;
Profile.epsilon=zeros(nbscan,2).*nan;
Profile.sh_fc=zeros(nbscan,2).*nan;
Profile.chi=zeros(nbscan,2).*nan;
Profile.tg_fc=zeros(nbscan,2).*nan;
Profile.tg_flag=zeros(nbscan,2).*nan;

% Profile.w=zeros(nbscan,1).*nan;
Profile.t=zeros(nbscan,1).*nan;
Profile.w=zeros(nbscan,1).*nan;
Profile.s=zeros(nbscan,1).*nan;
Profile.dnum=zeros(nbscan,1).*nan;

% loop along the pressure axis.
average_scan=@(x,y) (nanmean(x(y)));
for p=1:nbscan % p is the scan index.
    [~,indP] = sort(abs(Profile.P-Pr(p)));
    indP=indP(1);
    ind_ctdscan = indP-N_ctd/2:indP+N_ctd/2; % ind_scan is even
    scan.w   = average_scan(Profile.dPdt,ind_ctdscan(ind_ctdscan>0 & ind_ctdscan<LCTD));
    % check if the scan is not too shallow or too close to the end of the
    % profile. Also check if the speed if >20 cm s^{-1}
    if (ind_ctdscan(1)>0 && ind_ctdscan(end)<LCTD && scan.w>limit_speed) 
        ind_Pr_epsi = find(EpsiProfile.epsitime<CTDProfile.ctdtime(indP),1,'last');
        ind_epsiscan = ind_Pr_epsi-N_epsi/2:ind_Pr_epsi+N_epsi/2; % ind_scan is even

        % compute mean values of w,T,S of each scans
        scan.w    = average_scan(Profile.dPdt,ind_ctdscan);
        scan.t    = average_scan(CTDProfile.T,ind_ctdscan);
        scan.s    = average_scan(CTDProfile.S,ind_ctdscan);
        scan.pr   = average_scan(CTDProfile.P,ind_ctdscan); % I know this redondant with Pr axis
        scan.dnum = average_scan(CTDProfile.ctdtime,ind_ctdscan);
        scan.kvis=nu(scan.s,scan.t,scan.pr);
        scan.ktemp=kt(scan.s,scan.t,scan.pr);
        scan.kmax=fpump./scan.w;
        % compute spectra for acceleration channels.
        for c=[inda1 inda2 inda3]
            wh_channels=channels{c};
            scan.(wh_channels)=EpsiProfile.(wh_channels)(ind_epsiscan)*G; % time series in m.s^{-2}
            [scan.P.(wh_channels),fe] = pwelch(detrend(scan.(wh_channels)),nfft,[],nfft,df_epsi,'psd');
        end
        % get the filter transfer functions.
        if ~exist('h_freq','var')
            h_freq=get_filters_MADRE(Meta_Data,fe);
        end

        for c=1:length(channels)
            wh_channel=channels{c};
            switch wh_channel
                case 't1'
                    scan.(wh_channel)=EpsiProfile.(wh_channel)(ind_epsiscan).*dTdV(1); % time series in Celsius
                    [~,~,~,scan.chi(1),scan.tg_fc(1),scan.tg_flag(1)]=mod_efe_scan_chi(scan,wh_channel,Meta_Data,h_freq,FPO7noise);
                case 't2'
                    scan.(wh_channel)=EpsiProfile.(wh_channel)(ind_epsiscan).*dTdV(2); % time series in Celsius
                    [~,~,~,scan.chi(2),scan.tg_fc(2),scan.tg_flag(2)]=mod_efe_scan_chi(scan,wh_channel,Meta_Data,h_freq,FPO7noise);
                case 's1'
                    scan.(wh_channel)=EpsiProfile.(wh_channel)(ind_epsiscan).*twoG./(Sv1.*scan.w); % time series in m.s^{-1}
                    [~,~,~,~,~,scan.epsilon(1),scan.sh_fc(1),~]=mod_efe_scan_epsilon(scan,wh_channel,'a3',Meta_Data,h_freq);
                case 's2'
                    scan.(wh_channel)=EpsiProfile.(wh_channel)(ind_epsiscan).*twoG ./(Sv2.*scan.w); % time series in m.s^{-1}
                    [~,~,~,~,~,scan.epsilon(2),scan.sh_fc(2),~]=mod_efe_scan_epsilon(scan,wh_channel,'a3',Meta_Data,h_freq);
            end
        end
        Profile.epsilon(p,1)=scan.epsilon(1);
        Profile.epsilon(p,2)=scan.epsilon(2);
        Profile.sh_fc(p,1)=scan.sh_fc(1);
        Profile.sh_fc(p,2)=scan.sh_fc(2);

        Profile.chi(p,1)=scan.chi(1);
        Profile.chi(p,2)=scan.chi(2);
        Profile.tg_fc(p,1)=scan.tg_fc(1);
        Profile.tg_fc(p,2)=scan.tg_fc(2);
        Profile.tg_flag(p,1)=scan.tg_flag(1);
        Profile.tg_flag(p,2)=scan.tg_flag(2);

        Profile.w(p)=scan.w;
        Profile.t(p)=scan.t;
        Profile.s(p)=scan.s;
        Profile.dnum(p)=scan.dnum;
        
    end
end
Profile.fe=fe;

