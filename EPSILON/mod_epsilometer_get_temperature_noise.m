function [Tnoise1,Tnoise2]=mod_epsilometer_get_temperature_noise(Meta_Data,fc)
% [Tnoise1,Tnoise2]=mod_epsilometer_get_temperature_noise(Meta_Data,fc)
%
% adjust the epsilometer temperature channel electrical noise to the data by fitting the high freqeuncy part of the average T spectrum
% to the testbench electrical. It is used later to define the cut-off
% frequency for the TG spectrum integration.
%
%
% written by Arnaud Le Boyer 02/06/2020.
% edited by Nicole Couto 06/24/2020
%   - edited to use EpsiProfiles as they are just after raw data are
%   separated into profiles. This ensures that the whole deployment will be
%   used even if not all profiles are processed in L1 folder.


% Download EPSI and CTD profile
try
    load(fullfile(Meta_Data.paths.profiles,['Profiles_' Meta_Data.deployment '.mat']),'CTDProfiles','EpsiProfiles');
catch
    load(fullfile(Meta_Data.paths.profiles,['Profiles_' Meta_Data.deployment '.mat']),'CTDProfile','EpsiProfile');
    CTDProfiles = CTDProfile;
    EpsiProfiles = EpsiProfile;
end
switch Meta_Data.vehicle
    case 'FISH'
        CTD_Profiles = CTDProfiles.datadown;
        EPSI_Profiles = EpsiProfiles.datadown;
    case 'WW'
        CTD_Profiles = CTDProfiles.dataup;
        EPSI_Profiles = EpsiProfiles.dataup;
end


% Collect some additional data
fe=Meta_Data.PROCESS.fe;

% get FPO7 channel average noise to compute chi
switch Meta_Data.MAP.temperature
    case 'Tdiff'
        FPO7noise=load(fullfile(Meta_Data.paths.calibration,'FPO7_noise.mat'),'n0','n1','n2','n3');
    otherwise
        FPO7noise=load(fullfile(Meta_Data.paths.calibration,'FPO7_notdiffnoise.mat'),'n0','n1','n2','n3');
end

logf=log10(fe);
n0=FPO7noise.n0; n1=FPO7noise.n1; n2=FPO7noise.n2; n3=FPO7noise.n3;
noise=10.^(n0+n1.*logf+n2.*logf.^2+n3.*logf.^3);

nfft=Meta_Data.PROCESS.nfft;
Fs=Meta_Data.PROCESS.Fs_epsi;
dTdV(1)=Meta_Data.epsi.t1.dTdV; % define in mod_epsi_temperature_spectra
dTdV(2)=Meta_Data.epsi.t2.dTdV; % define in mod_epsi_temperature_spectra
dz = Meta_Data.PROCESS.dz;
limit_speed = 0.2;
tscan       = Meta_Data.PROCESS.tscan;
Fs_epsi     = Meta_Data.PROCESS.Fs_epsi;
N_epsi      = tscan.*Fs_epsi-mod(tscan*Fs_epsi,2);
Fs_ctd      = Meta_Data.PROCESS.Fs_ctd;
N_ctd       = tscan.*Fs_ctd-mod(tscan*Fs_ctd,2);

% Initialize variables
count1=1;
scan1=[];
allPt1=0;
allPt2=0;

for p=1:length(CTD_Profiles)
    
    
    % some cleaning on the Pressure channel
    % w = dPdt needs to be smooth.
    CTD_Profiles{p}.P = filloutliers(CTD_Profiles{p}.P,'center','movmedian',1000);
    
    CTD_Profiles{p} = structfun(@(x) fillmissing(x,'linear'),CTD_Profiles{p},'Un',0);
    EPSI_Profiles{p} = structfun(@(x) fillmissing(double(x),'linear'),EPSI_Profiles{p},'Un',0);
    
    %in case there is a mismatch between ctd and epsi time which still
    %happens as of April 17th 2020.
    EPSI_Profiles{p}.epsitime = EPSI_Profiles{p}.epsitime+ ...
        (CTD_Profiles{p}.ctdtime(1)-EPSI_Profiles{p}.epsitime(1));
    
    Prmin = min(CTD_Profiles{p}.P);
    Prmax = max(CTD_Profiles{p}.P);
    Profile = mod_epsilometer_merge_profile(CTD_Profiles{p},EPSI_Profiles{p},Prmin,Prmax);
    Pr = ceil(min(Profile.P)):dz:floor(max(Profile.P));
    nbscan = length(Pr);
    LCTD = length(Profile.P);
    
    
    for i=1:nbscan
        
        [~,indP]    = sort(abs(Profile.P-Pr(i)));
        indP        = indP(1);
        ind_ctdscan = indP-N_ctd/2:indP+N_ctd/2; % ind_scan is even
        
        scan1.w      = nanmean(Profile.dPdt(ind_ctdscan(ind_ctdscan>0 & ind_ctdscan<LCTD)));
        
        if ind_ctdscan(1)>N_ctd/2 && ind_ctdscan(end)<LCTD-N_ctd/2 && scan1.w>limit_speed
            
            scan1.Pr     = nanmean(Profile.P(ind_ctdscan));
            scan1.t      = nanmean(Profile.T(ind_ctdscan));
            scan1.s      = nanmean(Profile.S(ind_ctdscan));
            
            scan1.ktemp=kt(scan1.s,scan1.t,scan1.Pr);
            scan1 =mod_epsilometer_make_scan_v2(Profile,scan1,Meta_Data);
            % first get the spectrum in Volt so we can estimate the noise level and get
            % a cut-off freqeuncy
            [P1,~] = pwelch(detrend(scan1.t1)./dTdV(1),nfft,[],nfft,Fs,'psd');
            [P2,~] = pwelch(detrend(scan1.t2)./dTdV(1),nfft,[],nfft,Fs,'psd');
            
            allPt1=allPt1+P1;
            allPt2=allPt2+P2;
            count1=count1+1;
        end
    end %end loop through scans
    
end %end loop through profiles

noise1=allPt1./count1;
noise2=allPt2./count1;

coef1=nanmean(noise1(fe>fc)./noise(fe>fc));
coef2=nanmean(noise2(fe>fc)./noise(fe>fc));

polynoise1=polyfit(log10(fe(~isnan(noise))),log10(coef1.*noise(~isnan(noise))),3);

Tnoise1.n0=polynoise1(4);
Tnoise1.n1=polynoise1(3);
Tnoise1.n2=polynoise1(2);
Tnoise1.n3=polynoise1(1);

polynoise2=polyfit(log10(fe(~isnan(noise))),log10(coef2.*noise(~isnan(noise))),3);

Tnoise2.n0=polynoise2(4);
Tnoise2.n1=polynoise2(3);
Tnoise2.n2=polynoise2(2);
Tnoise2.n3=polynoise2(1);





