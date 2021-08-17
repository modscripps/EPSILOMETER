function Meta_Data=mod_epsi_temperature_spectra(Meta_Data,EPSI_Profile,CTD_Profile,title1,np,dsp,tscan)
%function Meta_Data=mod_epsi_temperature_spectra(Meta_Data,EPSI_Profile,CTD_Profile,title1,np,dsp,tscan,ctd_df)

% ALB Feb 2019.
%
% input
% Meta_Data: All environment informations
% EPSI and CTD data form 1 profile defined by EPSI_create_profiles 
% title = title of the plot
% np = Profile number. TODO: Could be directly implied when calling EPSI
% and CTD data
% tscan = length in seconds of the segment used to compute the spectra
% 
% ctd_df= sampling frequancy of the CTD instrument 
% ctd_df=6; % RBR ctd_df=16; % SBE

% default epsi_freq
epsi_df=Meta_Data.AFE.FS;

% define parameters to compute the spectra.
epsi_Lscan  = tscan*epsi_df;  
epsi_T      = length(EPSI_Profile.epsitime);
ctd_T       = length(CTD_Profile.ctdtime);

% make sure all the fields in the CTD structure are columns.
CTD_Profile=structfun(@(x) x(:),CTD_Profile,'un',0);
%Length_profile=86400*(EPSI_Profile.epsitime(end)-EPSI_Profile.epsitime(1));
Length_profile=EPSI_Profile.epsitime(end)-EPSI_Profile.epsitime(1); %NC - epsitime is now in seconds
ctd_df   = floor((ctd_T./Length_profile));
ctd_df   = Meta_Data.PROCESS.Fs_ctd;
ctd_Lscan  = tscan*ctd_df;  


fprintf('epsi time series %3.2f seconds.\n',epsi_T/epsi_df)
fprintf('ctd time series %3.2f seconds.\n',ctd_T/ctd_df)

% define number of segments.
nbscan  = floor(epsi_T/epsi_Lscan);

% define frequnecy axis for ctd and epsi spectra.
epsi_k=make_kaxis(tscan,epsi_df);
ctd_k=make_kaxis(tscan,ctd_df);

% compute fall rate
CTD_Profile.w = smoothdata([diff(CTD_Profile.P(:))./diff(CTD_Profile.ctdtime(:)*86400) ;nan],...
                           'movmean',10);

% we compute spectra on scan with 50% overlap.
nbscan=2*nbscan-1;
epsi_indscan = arrayfun(@(x) (1+floor(epsi_Lscan/2)*(x-1):1+floor(epsi_Lscan/2)*(x-1)+epsi_Lscan-1),1:nbscan,'un',0);
ctd_indscan = arrayfun(@(x) (1+floor(ctd_Lscan/2)*(x-1):1+floor(ctd_Lscan/2)*(x-1)+ctd_Lscan-1),1:nbscan,'un',0);
clear data_CTD

% split ctd data in segments. 
data_CTD = cell2mat(cellfun(@(x) CTD_Profile.T(x),ctd_indscan,'un',0)).';
% compute spectra for ctd data.
P11_ctd=alb_power_spectrum(data_CTD,1./tscan);

% initialize epsi data matrice
P11_epsi=zeros(2,nbscan,int32(epsi_Lscan));
% check if t1 and t2 channels exist. 
indt1=find(cellfun(@(x) strcmp(x,'t1'),Meta_Data.PROCESS.channels));
indt2=find(cellfun(@(x) strcmp(x,'t2'),Meta_Data.PROCESS.channels));
% define epsi data matrice after filling nans and removing outliers.
% compute the EPSI spectra
if ~isempty(indt1)
    EPSI_Profile.t1_volt=fillmissing(EPSI_Profile.t1_volt,'linear');
    EPSI_Profile.t1_volt=filloutliers(EPSI_Profile.t1_volt,'center','movmedian',1000);
    
    data_EPSI(indt1,:,:) = cell2mat(cellfun(@(x) filloutliers( ...
                EPSI_Profile.t1_volt(x),'center','movmedian',5), ...
                 epsi_indscan,'un',0)).';
    P11_epsi(indt1,:,:)=alb_power_spectrum(squeeze(data_EPSI(indt1,:,:)),1./tscan);
end
if ~isempty(indt2)
    EPSI_Profile.t2_volt=fillmissing(EPSI_Profile.t2_volt,'linear');
    EPSI_Profile.t2_volt=filloutliers(EPSI_Profile.t2_volt,'center','movmedian',1000);

    data_EPSI(indt2,:,:) = cell2mat(cellfun(@(x) filloutliers( ...
        EPSI_Profile.t2_volt(x),'center','movmedian',5),...
              epsi_indscan,'un',0)).';
    P11_epsi(indt2,:,:)=alb_power_spectrum(squeeze(data_EPSI(indt2,:,:)),1./tscan);
end

% split and median fall rate per segments. Used for the physical answer of
% the FPO7 probe (h_freq.FPO7).
w=cellfun(@(x) abs(median(CTD_Profile.w(x))),ctd_indscan);


% 1 side of the ctd spectra
ctd_indk=find(ctd_k>=0);
ctd_indk=ctd_indk(1:end-1);
ctd_k=ctd_k(ctd_indk);

% 1 side of the epsi spectra
epsi_indk=find(epsi_k>=0);
epsi_indk=epsi_indk(1:end-1);
epsi_k=epsi_k(epsi_indk);

% becasue we pick only 1 side I multiply by 2 
P11_ctd  = 2*squeeze(P11_ctd(:,ctd_indk));
% final ctd spectrum
P11_Tctd = nanmean(P11_ctd,1);


% becasue we pick only 1 side I multiply by 2 
P11_epsi = 2*squeeze(P11_epsi(:,:,epsi_indk));

% get transfer EPSI FPO7 transfert functions 
h_freq=get_filters_SOM(Meta_Data,epsi_k(:));
% compute fpo7 filters (they are speed dependent)
TFtemp=cell2mat(cellfun(@(x) h_freq.FPO7(x),num2cell(w),'un',0)).';

% apply the transfer function to correct the spectra 
% the ratio of the final ctd spectrum by the EPSI spectra 
% gives us the dTdV coeficient that allows to convert Volts to Celsius  
indsub1Hz=find(ctd_k<1); 
if ~isempty(indt1)
    P11_TFepsi(indt1,:,:)=squeeze(P11_epsi(indt1,:,:))./TFtemp;
    P11_T(indt1,:) = squeeze(nanmean(P11_TFepsi(indt1,:,:),2)); % Temperature gradient frequency spectra should be ?C^2/s^-2 Hz^-1 ????
    tempo_P11=interp1(epsi_k,P11_T(indt1,:),ctd_k);
    dTdV(1)=sqrt(nanmedian(P11_Tctd(indsub1Hz)./tempo_P11(indsub1Hz)));
    P11_T(indt1,:)= P11_T(indt1,:).*dTdV(1).^2;
end
if ~isempty(indt2)
    P11_TFepsi(indt2,:,:)=squeeze(P11_epsi(indt2,:,:))./TFtemp;
    P11_T(indt2,:) = squeeze(nanmean(P11_TFepsi(indt2,:,:),2)); % Temperature gradient frequency spectra should be ?C^2/s^-2 Hz^-1 ????
    tempo_P11=interp1(epsi_k,P11_T(indt2,:),ctd_k);
    dTdV(2)=sqrt(nanmedian(P11_Tctd(indsub1Hz)./tempo_P11(indsub1Hz)));
    P11_T(indt2,:)= P11_T(indt2,:).*dTdV(2).^2;
end


% A and B are the intermediary spectra. Only for plotting.
%A=squeeze(nanmean(P11_epsi(2,:,:),2)).*dTdV(1).^2;
A=squeeze(nanmean(P11_epsi(1,:,:),2));
AA=squeeze(nanmean(P11_epsi(1,:,:),2))./TFtemp.';
B=squeeze(nanmean(P11_epsi(1,:,:),2)).*dTdV(1).^2./h_freq.electFPO7.'.^2;

% Sensitivity of probe, nominal. Save dTdV in Meta Data it will be use in
% the batchprocess
Meta_Data.AFE.t1.cal=dTdV(1); 
Meta_Data.AFE.t2.cal=dTdV(2); 
save(fullfile(Meta_Data.datapath,'Meta_data.mat'),'Meta_Data');

% only for plotting: we getting the board noise
switch Meta_Data.AFE.temp_circuit
    case 'Tdiff'
        FPO7noise=load(fullfile(Meta_Data.CALIpath,'FPO7_noise.mat'),'n0','n1','n2','n3');
    otherwise
        FPO7noise=load(fullfile(Meta_Data.CALIpath,'FPO7_notdiffnoise.mat'),'n0','n1','n2','n3');
end
n0=FPO7noise.n0; n1=FPO7noise.n1; n2=FPO7noise.n2; n3=FPO7noise.n3;
logf=log10(1/3:1/3:160);
noise=10.^(n0+n1.*logf+n2.*logf.^2+n3.*logf.^3).*dTdV(1).^2;


% plot the spectra
close all
hold on
if ~isempty(indt1)
    loglog(epsi_k,P11_T(indt1,:),'c','linewidth',2)
else
    loglog(epsi_k,P11_T(indt2,:),'c','linewidth',2)
end
loglog(ctd_k,P11_Tctd,'r','linewidth',2)
%loglog(1/3:1/3:160,noise,'m','linewidth',2)
%loglog(epsi_k,10.^(mmp_noise).*dTdV(1).^2,'m--','linewidth',2)
loglog(epsi_k,A,'Color',.8* [1 1 1],'linewidth',2)
loglog(epsi_k,AA,'Color',.6* [1 1 1],'linewidth',2)
%loglog(epsi_k,B,'Color',.4* [1 1 1],'linewidth',2)
if ~isempty(indt1)
%    loglog(epsi_k,P11_T(indt1,:),'c','linewidth',2)
end
if ~isempty(indt2)
%    loglog(epsi_k,P11_T(indt2,:),'Color',.1* [1 1 1],'linewidth',2)
end

set(gca,'XScale','log','YScale','log')
xlabel('Hz','fontsize',20)
ylabel('C^2/Hz','fontsize',20)
title1=sprintf('%s cast %i - temperature',title1,np);
title(title1,'fontsize',20)
%legend('CTD','noise','raw','t1./TF_{elec}','t1','t2','location','southwest')
legend('t1','CTD','raw(Volt^2/Hz)','corrected (Volt^2/Hz)','location','northeast')
grid on
ylim([1e-13 1])
xlim([1/15 170])
set(gca,'fontsize',20)
fig=gcf;fig.PaperPosition=[0 0 8 6];
filename=sprintf('%sTctd_Tepsi_comp_cast%i_t%i.png',Meta_Data.L1path,np,1);
figureStamp(mfilename('fullpath'))
print('-dpng2',filename)
if dsp==1
    ind_t=1;
    for i=1:nbscan
        close all
        A=squeeze(P11_epsi(ind_t,i,:)).*dTdV(ind_t).^2;
        B=squeeze(P11_epsi(ind_t,i,:)).*dTdV(ind_t).^2./h_freq.electFPO7.'.^2;
        loglog(epsi_k,A,'c')
        hold on
        loglog(epsi_k,10.^(noise).*dTdV(ind_t).^2,'m','linewidth',2)
        loglog(epsi_k,10.^(mmp_noise).*dTdV(ind_t).^2,'m--','linewidth',2)
        %loglog(epsi_k,A,'Color',.6* [1 1 1],'linewidth',2)
        %loglog(epsi_k,B,'Color',.4* [1 1 1],'linewidth',2)
        loglog(epsi_k,squeeze(P11_TFepsi(1,i,:)),'r','linewidth',2)
        loglog(ctd_k,P11_ctd(i,:),'k')
        set(gca,'XScale','log','YScale','log')
        legend('t1 raw','noise',' mmp noise','t1','SBE')
        xlabel('Hz','fontsize',20)
        ylabel('C^2/Hz','fontsize',20)
        title(title1,'fontsize',20)
        set(gca,'fontsize',20)
        ylim([1e-13 1])
        xlim([1/15 170])
        grid on
        pause
    end
close all
end


