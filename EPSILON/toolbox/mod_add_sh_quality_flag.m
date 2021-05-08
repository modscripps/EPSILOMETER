function Profile=mod_add_sh_quality_flag(Profile,Meta_Data,H1,H2,fH,fc1,fc2)

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
df=Meta_Data.df_epsi;

tscan=Meta_Data.tscan;

%% Gravity  ... of the situation :)
G       = 9.81;
twoG= 2*G;
%% define a Pressure axis to an which I will compute epsilon and chi. 
%  The spectra will be nfft long centered around P(z) +/- tscan/2. 
%  

Pr=Profile.pr;
nbscan=length(Pr);

% number of samples for a scan. I make sure it is always even
N_pr=tscan.*df-mod(tscan*df,2);

%initialize process flags
Profile.sh_qcflag=zeros(nbscan,2).*nan;
TFnoise=@(x,y) (interp1(fH,y,x));

for p=1:nbscan % p is the scan index.
    [~,indP] = sort(abs(Profile.P-Pr(p)));
    indP=indP(1);
    ind_scan = indP-N_pr/2:indP+N_pr/2; % ind_scan is even
    
    wh_channels=channels{inda3};
    scan.(wh_channels)=Profile.(wh_channels)(ind_scan)*G; % time series in m.s^{-2}
    [scan.P.(wh_channels),fe] = pwelch(detrend(scan.(wh_channels)),nfft,[],nfft,df_epsi,'psd');
    
    u1_vibration=TFnoise(scan.P.(wh_channels),H1);
    u2_vibration=TFnoise(scan.P.(wh_channels),H2);

    for c=1:nb_channels
        wh_channel=channels{c};
        switch wh_channel
            case 's1'
                scan.(wh_channel)=Profile.(wh_channel)(ind_scan).*twoG./(Sv1.*scan.w); % time series in m.s^{-1}
                [P1,~,~,~,~,~,~,fe]=mod_efe_scan_epsilon(scan,wh_channel,'a3',Meta_Data,h_freq);

            case 's2'
                scan.(wh_channel)=Profile.(wh_channel)(ind_scan).*twoG ./(Sv2.*scan.w); % time series in m.s^{-1}
                [P2,~,~,~,~,~,~,fe]=mod_efe_scan_epsilon(scan,wh_channel,'a3',Meta_Data,h_freq);
        end
    end
    
    % get the ratio of U_obs / U_vibration (from low epsilon region)
    qc_flag1= log10(P1./ smoothdata(u1_vibration,'movmean',5));
    qc_flag2= log10(P2./ smoothdata(u2_vibration,'movmean',5));

    qc_flag1=nanmean(qc_flag1(fe>fc1 & fe<fc2));
    qc_flag2=nanmean(qc_flag2(fe>fc1 & fe<fc2));

    Profile.sh_qcflag(p,:)=[qc_flag1 qc_flag2];
end
