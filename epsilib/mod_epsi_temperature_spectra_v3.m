function Meta_Data=mod_epsi_temperature_spectra_v3(Meta_Data,Profile)
%function Meta_Data=mod_epsi_temperature_spectra_v3(Meta_Data,Profile,Profile,title1,np,dsp,tscan)
%
% NC edited May 2021 to take combined Profile structure as input, instead
% of separate EPSI_Profile and CTD_Profile.
%   - PROS: This allows user to calibrate CTD using just one Profile or
%   timeseries without having to process the whole deployment
%   - CONS: The chosen profile might not be the best chunk of data for the
%   calibration
%  - SOLUTION: Find the longest profile first, before you call this. This
%  is done in epsi_class f_calibrateTemperature
%  - TODO: If only Meta_Data is input, find the longest profile within this
%  function and do the calibration on that
%
% ALB Feb 2019.
%
% input
% Meta_Data: All environment informations
% EPSI and CTD data form 1 profile defined by EPSI_create_profiles
% tscan = length in seconds of the segment used to compute the spectra
%
% ctd_df= sampling frequancy of the CTD instrument
% ctd_df=6; % RBR ctd_df=16; % SBE

% For the temperature correction, we want a good chunk of data. Try
% tscan = 50 or if the profile is too short, find the longest profile
% in the deployment and set tscan to about 0.5-0.8 times the profile
% length

% default epsi and ctd sampling frequencies
try
    epsi_df=Meta_Data.AFE.FS;
catch
    epsi_df=Meta_Data.PROCESS.Fs_epsi; %Check this,  I think it's wrong
end

ctd_df = Meta_Data.PROCESS.Fs_ctd;

% id profile
tscanAlternate = floor(0.8*length(Profile.epsi.time_s)/epsi_df);
tscanDefault = 20;
tscan = min([tscanDefault,tscanAlternate]);

% define parameters to compute the spectra.
epsi_Lscan  = tscan*epsi_df;
epsi_T      = length(Profile.epsi.time_s);
ctd_Lscan   = tscan*ctd_df;
ctd_T       = length(Profile.ctd.time_s);


% make sure all the fields in the CTD structure are columns.
Profile=structfun(@(x) x(:),Profile,'un',0);
%Length_profile=86400*(Profile.epsi.time_s(end)-Profile.epsi.time_s(1));
Length_profile=Profile.epsi.time_s(end)-Profile.epsi.time_s(1); %NC - epsitime is now in seconds
% ctd_df   = floor((ctd_T./Length_profile));
% ctd_df   = Meta_Data.PROCESS.Fs_ctd;
ctd_Lscan  = tscan*ctd_df;

fprintf('epsi time series %3.2f seconds.\n',epsi_T/epsi_df)
fprintf('ctd time series %3.2f seconds.\n',ctd_T/ctd_df)

% define number of segments.
% NC - Added nbscan for ctd and compare the two. Sometimes it was making
% one two many scans for ctd to use based on epsi
nbscan_epsi  = floor(epsi_T/epsi_Lscan);
nbscan_ctd   = floor(ctd_T/ctd_Lscan);
nbscan = min([nbscan_epsi,nbscan_ctd]);

% define frequnecy axis for ctd and epsi spectra.
epsi_k=make_kaxis(tscan,epsi_df);
ctd_k=make_kaxis(tscan,ctd_df);

% compute fall rate
% Profile.w = smoothdata([diff(Profile.z(:))./diff(Profile.ctd.time_s(:)*86400) ;nan],...
%                            'movmean',10
% ctdtime is in seconds
Profile.w = smoothdata([diff(Profile.ctd.z(:))./diff(Profile.ctd.time_s(:)) ;nan],...
                           'movmean',10);

% we compute spectra on scan with 50% overlap.
nbscan=2*nbscan-1;
epsi_indscan = arrayfun(@(x) (1+floor(epsi_Lscan/2)*(x-1):1+floor(epsi_Lscan/2)*(x-1)+epsi_Lscan-1),1:nbscan,'un',0);
ctd_indscan = arrayfun(@(x) (1+floor(ctd_Lscan/2)*(x-1):1+floor(ctd_Lscan/2)*(x-1)+ctd_Lscan-1),1:nbscan,'un',0);
clear data_CTD

% split ctd data in segments.
data_CTD = cell2mat(cellfun(@(x) Profile.ctd.T(x),ctd_indscan,'un',0)).';
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
    Profile.epsi.t1_volt=fillmissing(Profile.epsi.t1_volt,'linear');
    Profile.epsi.t1_volt=filloutliers(Profile.epsi.t1_volt,'center','movmedian',1000);

    data_EPSI(indt1,:,:) = cell2mat(cellfun(@(x) filloutliers( ...
                Profile.epsi.t1_volt(x),'center','movmedian',5), ...
                 epsi_indscan,'un',0)).';
    P11_epsi(indt1,:,:)=alb_power_spectrum(squeeze(data_EPSI(indt1,:,:)),1./tscan);
end
if ~isempty(indt2)
    Profile.epsi.t2_volt=fillmissing(Profile.epsi.t2_volt,'linear');
    Profile.epsi.t2_volt=filloutliers(Profile.epsi.t2_volt,'center','movmedian',1000);

    data_EPSI(indt2,:,:) = cell2mat(cellfun(@(x) filloutliers( ...
        Profile.epsi.t2_volt(x),'center','movmedian',5),...
              epsi_indscan,'un',0)).';
    P11_epsi(indt2,:,:)=alb_power_spectrum(squeeze(data_EPSI(indt2,:,:)),1./tscan);
end

% split and median fall rate per segments. Used for the physical answer of
% the FPO7 probe (h_freq.FPO7).
w=cellfun(@(x) abs(median(Profile.ctd.dzdt(x))),ctd_indscan);

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
A1=squeeze(nanmean(P11_epsi(indt1,:,:),2));
AA1=squeeze(nanmean(P11_epsi(indt1,:,:),2))./TFtemp.';
B1=squeeze(nanmean(P11_epsi(indt1,:,:),2)).*dTdV(1).^2./h_freq.electFPO7.'.^2;

A2=squeeze(nanmean(P11_epsi(indt1,:,:),2));
AA2=squeeze(nanmean(P11_epsi(indt1,:,:),2))./TFtemp.';
B2=squeeze(nanmean(P11_epsi(indt1,:,:),2)).*dTdV(1).^2./h_freq.electFPO7.'.^2;

% Sensitivity of probe, nominal. Save dTdV in Meta Data it will be use in
% the batchprocess
Meta_Data.AFE.t1.cal=dTdV(1);
Meta_Data.AFE.t2.cal=dTdV(2);
save(fullfile(Meta_Data.paths.data,'Meta_data.mat'),'Meta_Data');

% only for plotting: we getting the board noise
try
switch Meta_Data.AFE.temp_circuit
    case 'Tdiff'
        FPO7noise=load(fullfile(Meta_Data.paths.calibration,'FPO7_noise.mat'),'n0','n1','n2','n3');
    otherwise
        FPO7noise=load(fullfile(Meta_Data.paths.calibration,'FPO7_notdiffnoise.mat'),'n0','n1','n2','n3');
end
catch
switch Meta_Data.MAP.temperature
    case 'Tdiff'
        FPO7noise=load(fullfile(Meta_Data.paths.calibration,'FPO7_noise.mat'),'n0','n1','n2','n3');
    otherwise
        FPO7noise=load(fullfile(Meta_Data.paths.calibration,'FPO7_notdiffnoise.mat'),'n0','n1','n2','n3');
end
end
n0=FPO7noise.n0; n1=FPO7noise.n1; n2=FPO7noise.n2; n3=FPO7noise.n3;
logf=log10(1/3:1/3:epsi_df/2);
noise=10.^(n0+n1.*logf+n2.*logf.^2+n3.*logf.^3).*dTdV(1).^2;

% plot the spectra
close all

% t1
% ---------------------------------
figure
hold on
if ~isempty(indt1)
    loglog(epsi_k,P11_T(indt1,:),'c','linewidth',2)
end
loglog(ctd_k,P11_Tctd,'r','linewidth',2)
%loglog(1/3:1/3:160,noise,'m','linewidth',2)
%loglog(epsi_k,10.^(mmp_noise).*dTdV(1).^2,'m--','linewidth',2)
loglog(epsi_k,A1,'Color',.8* [1 1 1],'linewidth',2)
loglog(epsi_k,AA1,'Color',.6* [1 1 1],'linewidth',2)
%loglog(epsi_k,B,'Color',.4* [1 1 1],'linewidth',2)
if ~isempty(indt1)
%    loglog(epsi_k,P11_T(indt1,:),'c','linewidth',2)
end

set(gca,'XScale','log','YScale','log')
xlabel('Hz','fontsize',20)
ylabel('C^2/Hz','fontsize',20)
% Title for figure
titleStr = strrep([Meta_Data.mission ' ' Meta_Data.vehicle_name ' ' Meta_Data.deployment],'_','\_');
if isfield(Profile,'profNum')
titleStr={sprintf('%s cast %i - temperature_ ',titleStr,Profile.profNum);...
                sprintf('dTdV = %2.2f',dTdV(1))};
else
  titleStr={['epsitime ' datestr(Profile.epsi.dnum(1),'HH:MM:SS') ' - ' datestr(Profile.epsi.dnum(end),'HH:MM:SS')]; ...
      sprintf('dTdV = %2.2f',dTdV(1))};
end
title(titleStr,'fontsize',20)
%legend('CTD','noise','raw','t1./TF_{elec}','t1','t2','location','southwest')
legend('t1','CTD','raw(Volt^2/Hz)','corrected (Volt^2/Hz)','location','northeast')
grid on
ylim([1e-13 1])
xlim([1/15 170])
set(gca,'fontsize',20)
fig=gcf;fig.PaperPosition=[0 0 8 6];
if isfield(Profile,'profNum')
    filename=fullfile(Meta_Data.paths.figures,sprintf('Tctd_Tepsi_comp_cast%i_t1.png',Profile.profNum));
else
    filename=fullfile(Meta_Data.paths.figures,'Tctd_Tepsi_comp_t1.png');
end
figureStamp(mfilename('fullpath'))
print('-dpng2',filename)


% t2
% ---------------------------------
figure
hold on
if ~isempty(indt2)
    loglog(epsi_k,P11_T(indt2,:),'c','linewidth',2)
end
loglog(ctd_k,P11_Tctd,'r','linewidth',2)
%loglog(1/3:1/3:160,noise,'m','linewidth',2)
%loglog(epsi_k,10.^(mmp_noise).*dTdV(1).^2,'m--','linewidth',2)
loglog(epsi_k,A2,'Color',.8* [1 1 1],'linewidth',2)
loglog(epsi_k,AA2,'Color',.6* [1 1 1],'linewidth',2)
%loglog(epsi_k,B,'Color',.4* [1 1 1],'linewidth',2)
if ~isempty(indt2)
%    loglog(epsi_k,P11_T(indt2,:),'Color',.1* [1 1 1],'linewidth',2)
end

set(gca,'XScale','log','YScale','log')
xlabel('Hz','fontsize',20)
ylabel('C^2/Hz','fontsize',20)
% Title for figure
titleStr = strrep([Meta_Data.mission ' ' Meta_Data.vehicle_name ' ' Meta_Data.deployment],'_','\_');
if isfield(Profile,'profNum')
titleStr={sprintf('%s cast %i - temperature_ ',titleStr,Profile.profNum);...
                sprintf('dTdV = %2.2f',dTdV(2))};
else
  titleStr={['epsitime ' datestr(Profile.epsi.dnum(1),'HH:MM:SS') ' - ' datestr(Profile.epsi.dnum(end),'HH:MM:SS')]; ...
      sprintf('dTdV = %2.2f',dTdV(2))};
end
title(titleStr,'fontsize',20)
%legend('CTD','noise','raw','t1./TF_{elec}','t1','t2','location','southwest')
legend('t2','CTD','raw(Volt^2/Hz)','corrected (Volt^2/Hz)','location','northeast')
grid on
ylim([1e-13 1])
xlim([1/15 170])
set(gca,'fontsize',20)
fig=gcf;fig.PaperPosition=[0 0 8 6];
if isfield(Profile,'profNum')
    filename=fullfile(Meta_Data.paths.figures,sprintf('Tctd_Tepsi_comp_cast%i_t2.png',Profile.profNum));
else
    filename=fullfile(Meta_Data.paths.figures,'Tctd_Tepsi_comp_t2.png');
end
figureStamp(mfilename('fullpath'))
print('-dpng2',filename)






% if dsp==1
%     ind_t=1;
%     for i=1:nbscan
%         close all
%         A=squeeze(P11_epsi(ind_t,i,:)).*dTdV(ind_t).^2;
%         B=squeeze(P11_epsi(ind_t,i,:)).*dTdV(ind_t).^2./h_freq.electFPO7.'.^2;
%         loglog(epsi_k,A,'c')
%         hold on
%         loglog(epsi_k,10.^(noise).*dTdV(ind_t).^2,'m','linewidth',2)
%         loglog(epsi_k,10.^(mmp_noise).*dTdV(ind_t).^2,'m--','linewidth',2)
%         %loglog(epsi_k,A,'Color',.6* [1 1 1],'linewidth',2)
%         %loglog(epsi_k,B,'Color',.4* [1 1 1],'linewidth',2)
%         loglog(epsi_k,squeeze(P11_TFepsi(1,i,:)),'r','linewidth',2)
%         loglog(ctd_k,P11_ctd(i,:),'k')
%         set(gca,'XScale','log','YScale','log')
%         legend('t1 raw','noise',' mmp noise','t1','SBE')
%         xlabel('Hz','fontsize',20)
%         ylabel('C^2/Hz','fontsize',20)
%         title(titleStr,'fontsize',20)
%         set(gca,'fontsize',20)
%         ylim([1e-13 1])
%         xlim([1/15 170])
%         grid on
%         pause
%     end
% close all
% end
