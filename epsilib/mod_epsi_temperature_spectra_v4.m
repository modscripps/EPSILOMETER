function Meta_Data=mod_epsi_temperature_spectra_v4(Meta_Data,Profile,saveData,plotData)
%function Meta_Data=mod_epsi_temperature_spectra_v3(Meta_Data,Profile,Profile,title1,np,dsp,tscan)
%
% ALB Sept 1st 2024 
% With the new version of the epsilib I need to save dTdV in a text files
% (similar format as Sv for the Sher probes). This is great becasue we will
% be able to keep track to the dTdV.... which should change a function of
% CTD so I need to save that too in the text file. 
%
% NC Sept 2021 - added dTdV_profile that gives dTdV for each scan length
% (~50 seconds)
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

% Options for saving and plotting data
if nargin<4
    plotData = 1;
    if nargin<3
        saveData = 1;
    end
end
    
% default epsi and ctd sampling frequencies
try
    epsi_df=Meta_Data.AFE.FS;
catch
    epsi_df=Meta_Data.PROCESS.Fs_epsi;
end

ctd_df = Meta_Data.PROCESS.Fs_ctd;

% id profile
tscanAlternate = floor(0.8*length(Profile.epsi.time_s)/epsi_df);
tscanDefault = 50;
tscan = min([tscanDefault,tscanAlternate/2]); %NC - alternate length is half the profile length, so you can have some overlap

% define parameters to compute the spectra.
epsi_Lscan  = tscan*epsi_df;
epsi_T      = length(Profile.epsi.time_s);
ctd_Lscan   = tscan*ctd_df;
ctd_T       = length(Profile.ctd.time_s);


% make sure all the fields in the CTD structure are columns.
% NC 9/13/21 - Comment out. Profile strucutre has changed since this line
% was written
%Profile=structfun(@(x) x(:),Profile,'un',0);
%Length_profile=86400*(Profile.epsi.time_s(end)-Profile.epsi.time_s(1));
epsi_profile_duration=Profile.epsi.time_s(end)-Profile.epsi.time_s(1); %NC - epsitime is now in seconds
ctd_profile_duration=Profile.ctd.time_s(end)-Profile.ctd.time_s(1); %NC - epsitime is now in seconds

% ctd_df   = floor((ctd_T./Length_profile));
% ctd_df   = Meta_Data.PROCESS.Fs_ctd;
% ctd_Lscan  = tscan*ctd_df;

% NC 9/29/21 - Comment out display of timeseries lengths
% fprintf('epsi time series %3.2f seconds.\n',epsi_profile_duration)
% fprintf('ctd time series %3.2f seconds.\n',ctd_profile_duration)

% define number of segments.
% NC - Added nbscan for ctd and compare the two. Sometimes it was making
% one too many scans for ctd to use based on epsi
nbscan_epsi  = floor(epsi_T/epsi_Lscan);
nbscan_ctd   = floor(ctd_T/ctd_Lscan);
nbscan = min([nbscan_epsi,nbscan_ctd]);

% % compute fall rate
% % Profile.w = smoothdata([diff(Profile.z(:))./diff(Profile.ctd.time_s(:)*86400) ;nan],...
% %                            'movmean',10
% % ctdtime is in seconds
% Profile.w = smoothdata([diff(Profile.ctd.z(:))./diff(Profile.ctd.time_s(:)) ;nan],...
%                            'movmean',10);

% we compute spectra on scan with 50% overlap.
if nbscan>1
nbscan=2*nbscan-1;
epsi_indscan = arrayfun(@(x) (1+floor(epsi_Lscan/2)*(x-1):1+floor(epsi_Lscan/2)*(x-1)+epsi_Lscan-1),1:nbscan,'un',0);
ctd_indscan = arrayfun(@(x) (1+floor(ctd_Lscan/2)*(x-1):1+floor(ctd_Lscan/2)*(x-1)+ctd_Lscan-1),1:nbscan,'un',0);
clear data_CTD
T=filloutliers(Profile.ctd.T,'linear','movmean',100);

% split ctd data in segments.
data_CTD = cell2mat(cellfun(@(x) T(x),ctd_indscan,'un',0)).';
time_CTD = cell2mat(cellfun(@(x) Profile.ctd.time_s(x),ctd_indscan,'un',0)).';
pr_CTD = cell2mat(cellfun(@(x) Profile.ctd.P(x),ctd_indscan,'un',0)).';

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
time_EPSI = cell2mat(cellfun(@(x) filloutliers( ...
                Profile.epsi.time_s(x),'center','movmedian',5), ...
                 epsi_indscan,'un',0)).';

% split and median fall rate per segments. Used for the physical answer of
% the FPO7 probe (h_freq.FPO7).
dzdt_smooth = movmean(Profile.ctd.dzdt,epsi_df);
w=cellfun(@(x) abs(median(dzdt_smooth(x))),ctd_indscan);

% define frequnecy axis for ctd and epsi spectra.
epsi_f=make_kaxis(tscan,epsi_df);
ctd_f=make_kaxis(tscan,ctd_df);

% 1 side of the ctd spectra
ctd_indf=find(ctd_f>=0);
ctd_indf=ctd_indf(1:end-1);
ctd_f=ctd_f(ctd_indf);

% 1 side of the epsi spectra
epsi_indf=find(epsi_f>=0);
epsi_indf=epsi_indf(1:end-1);
epsi_f=epsi_f(epsi_indf);

% becasue we pick only 1 side I multiply by 2
P11_ctd  = 2*squeeze(P11_ctd(:,ctd_indf));
% final ctd spectrum
P11_ctd_mean = nanmean(P11_ctd,1);

% becasue we pick only 1 side I multiply by 2
P11_epsi = 2*squeeze(P11_epsi(:,:,epsi_indf));

% get transfer EPSI FPO7 transfert functions
if isfield(Meta_Data,'MADRE')
    h_freq = get_filters_MADRE(Meta_Data,epsi_f(:));
else
h_freq=get_filters_SOM(Meta_Data,epsi_f(:));
end
% compute fpo7 filters (they are speed dependent)
TFtemp=cell2mat(cellfun(@(x) h_freq.FPO7(x),num2cell(w),'un',0)).';

% apply the transfer function to correct the spectra
% the ratio of the final ctd spectrum by the EPSI spectra
% gives us the dTdV coeficient that allows to convert Volts to Celsius
indsub1Hz = ctd_f<1;

% TODO 9/14/21 - This only works if there is more than one scan. Make it
% work if there is exactly one scan too. 
if nbscan>=1 

    if ~isempty(indt1)
        P11_epsi_TF(indt1,:,:)=squeeze(P11_epsi(indt1,:,:))./TFtemp;
        P11_epsi_TF_interp(indt1,:,:)=interp1(epsi_f,squeeze(P11_epsi_TF(indt1,:,:)).',ctd_f).';
        if ndims(P11_epsi_TF_interp)==3
            dTdV_profile(indt1,:) = sqrt(nanmedian(P11_ctd(:,indsub1Hz)./squeeze(P11_epsi_TF_interp(indt1,:,indsub1Hz)),2));
        elseif ndims(P11_epsi_TF_interp)==2
            dTdV_profile(indt1,:) = sqrt(nanmedian(P11_ctd(:,indsub1Hz)./P11_epsi_TF_interp(indt1,indsub1Hz)));
        end
        P11_epsi_TF_mean(indt1,:) = squeeze(nanmean(P11_epsi_TF(indt1,:,:),2)); % Temperature gradient frequency spectra should be ?C^2/s^-2 Hz^-1 ????
        P11_epsi_TF_mean_interp(indt1,:)=interp1(epsi_f,P11_epsi_TF_mean(indt1,:),ctd_f);
        dTdV(1)=sqrt(nanmedian(P11_ctd_mean(indsub1Hz)./P11_epsi_TF_mean_interp(indt1,indsub1Hz)));
        P11_epsi_TF_mean_calibrated(indt1,:)= P11_epsi_TF_mean(indt1,:).*dTdV(1).^2;
    end
    if ~isempty(indt2)
        P11_epsi_TF(indt2,:,:)=squeeze(P11_epsi(indt2,:,:))./TFtemp;
        P11_epsi_TF_interp(indt2,:,:)=interp1(epsi_f,squeeze(P11_epsi_TF(indt2,:,:)).',ctd_f).';
        if ndims(P11_epsi_TF_interp)==3
            dTdV_profile(indt2,:) = sqrt(nanmedian(P11_ctd(:,indsub1Hz)./squeeze(P11_epsi_TF_interp(indt2,:,indsub1Hz)),2));
        elseif ndims(P11_epsi_TF_interp)==2
            dTdV_profile(indt2,:) = sqrt(nanmedian(P11_ctd(:,indsub1Hz)./P11_epsi_TF_interp(indt2,indsub1Hz)));
        end
        P11_epsi_TF_mean(indt2,:) = squeeze(nanmean(P11_epsi_TF(indt2,:,:),2)); % Temperature gradient frequency spectra should be ?C^2/s^-2 Hz^-1 ????
        P11_epsi_TF_mean_interp(indt2,:)=interp1(epsi_f,P11_epsi_TF_mean(indt2,:),ctd_f);
        dTdV(2)=sqrt(nanmedian(P11_ctd_mean(indsub1Hz)./P11_epsi_TF_mean_interp(indt2,indsub1Hz)));
        P11_epsi_TF_mean_calibrated(indt2,:)= P11_epsi_TF_mean(indt2,:).*dTdV(2).^2;
    end

else 
    dTdV = [nan nan];
    dTdV_profile = nan(2,nbscan);
end


% A and B are the intermediary spectra. Only for plotting.
if ~isempty(indt1)
    A1=squeeze(nanmean(P11_epsi(indt1,:,:),2));
    AA1=squeeze(nanmean(P11_epsi(indt1,:,:),2))./TFtemp.';
    B1=squeeze(nanmean(P11_epsi(indt1,:,:),2)).*dTdV(1).^2./h_freq.electFPO7.'.^2;
end

if ~isempty(indt2)
    A2=squeeze(nanmean(P11_epsi(indt2,:,:),2));
    AA2=squeeze(nanmean(P11_epsi(indt2,:,:),2))./TFtemp.';
    B2=squeeze(nanmean(P11_epsi(indt2,:,:),2)).*dTdV(2).^2./h_freq.electFPO7.'.^2;
end
% Sensitivity of probe, nominal. Save dTdV in Meta Data it will be use in
% the batchprocess
if isfield(Meta_Data,'AFE')
    field_name = 'AFE';
elseif isfield(Meta_Data,'epsi')
    field_name = 'epsi';
end
Meta_Data.(field_name).t1.cal=dTdV(1);
Meta_Data.(field_name).t2.cal=dTdV(2);
% NC 9/14/21 - Added dTdV_profile and rangeT to try to come up with a
% better way to compute dTdV for a deployment
Meta_Data.(field_name).t1.cal_profile=dTdV_profile(1,:);
try
    Meta_Data.(field_name).t2.cal_profile=dTdV_profile(2,:);
catch
    Meta_Data.(field_name).t2.cal_profile=dTdV_profile(1,:);
end
Meta_Data.(field_name).t1.ctd_Tmin = nanmin(data_CTD.');
Meta_Data.(field_name).t1.ctd_Tmax = nanmax(data_CTD.');
Meta_Data.(field_name).t1.pr = nanmean(pr_CTD.');


% Save dTdV in Meta_Data
if saveData
    save(fullfile(Meta_Data.paths.data,'Meta_data.mat'),'Meta_Data');
    % TODO change this so we append the new dVdT.
    % Right now I am doing a quick an dirty open write save.
    %t1
    calibration_file =  ...
        sprintf('Calibration_%s.txt',Meta_Data.(field_name).t1.SN);
    calibration_path =  ...
        fullfile(Meta_Data.AFE.tempcal_path,Meta_Data.(field_name).t1.SN);

    if ~exist(calibration_path,'dir')
       mkdir(calibration_path);
    end
    fid=fopen(fullfile(calibration_path,calibration_file),"w");
    calibration_string=sprintf('%s\r\n %03s, %2.1f, %s\r\n',...
                        datestr(now),Meta_Data.(field_name).t1.SN, ...
                        dTdV(1), ...
                        Meta_Data.CTD.cal.SN);
    fwrite(fid,calibration_string);
    fclose(fid);
    % t2
    calibration_file =  ...
        sprintf('Calibration_%s.txt',Meta_Data.(field_name).t2.SN);
    calibration_path =  ...
        fullfile(Meta_Data.AFE.tempcal_path,Meta_Data.(field_name).t2.SN);

    if ~exist(calibration_path,'dir')
       mkdir(calibration_path);
    end
    fid=fopen(fullfile(calibration_path,calibration_file),"w");
    calibration_string=sprintf('%s\r\n %03s, %2.1f, %s\r\n',...
                        datestr(now),Meta_Data.(field_name).t2.SN, ...
                        dTdV(2), ...
                        Meta_Data.CTD.cal.SN);
    fwrite(fid,calibration_string);
    fclose(fid);
    
end

% Plot data
if  plotData
    % only for plotting: we getting the board noise
    switch Meta_Data.(field_name).temp_circuit
        case 'Tdiff'
            FPO7noise=load(fullfile(Meta_Data.paths.calibration,'FPO7_noise.mat'),'n0','n1','n2','n3');
        otherwise
            FPO7noise=load(fullfile(Meta_Data.paths.calibration,'FPO7_notdiffnoise.mat'),'n0','n1','n2','n3');
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
        p1 = loglog(epsi_f,P11_epsi_TF_mean_calibrated(indt1,:),'c','linewidth',2,'displayname','t1');
    end
    p2 = loglog(ctd_f,P11_ctd_mean,'r','linewidth',2,'displayname','CTD');
    %loglog(1/3:1/3:160,noise,'m','linewidth',2)
    %loglog(epsi_k,10.^(mmp_noise).*dTdV(1).^2,'m--','linewidth',2)
    try
    p3 = loglog(epsi_f,A1,'Color',.8* [1 1 1],'linewidth',2,'displayname','raw(Volt^2/Hz)');
    p4 = loglog(epsi_f,AA1,'Color',.6* [1 1 1],'linewidth',2,'displayname','corrected (Volt^2/Hz)');
    catch
    p3 = loglog(epsi_f,A2,'Color',.8* [1 1 1],'linewidth',2,'displayname','raw(Volt^2/Hz)');
    p4 = loglog(epsi_f,AA2,'Color',.6* [1 1 1],'linewidth',2,'displayname','corrected (Volt^2/Hz)');
    end
    %loglog(epsi_k,B,'Color',.4* [1 1 1],'linewidth',2)
    if ~isempty(indt1)
        %    loglog(epsi_k,P11_T(indt1,:),'c','linewidth',2)
    end
    
    set(gca,'XScale','log','YScale','log')
    xlabel('Hz','fontsize',20)
    ylabel('C^2/Hz','fontsize',20)

    % Title for figure
    titleStr = strrep([Meta_Data.cruise_name ' ' Meta_Data.vehicle_name ' ' Meta_Data.deployment],'_','\_');
    if isfield(Profile,'profNum')
        titleStr={sprintf('%s cast %i - temperature_ ',titleStr,Profile.profNum);...
            sprintf('dTdV = %2.2f',dTdV(1))};
    else
        titleStr={['epsitime ' datestr(Profile.epsi.dnum(1),'HH:MM:SS') ' - ' datestr(Profile.epsi.dnum(end),'HH:MM:SS')]; ...
            sprintf('dTdV = %2.2f',dTdV(1))};
    end
    title(titleStr,'fontsize',20)

    legend([p1,p2,p3,p4(1)],'location','northeast')
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
        p1 = loglog(epsi_f,P11_epsi_TF_mean_calibrated(indt2,:),'c','linewidth',2,'displayname','t2');
    end
    p2 = loglog(ctd_f,P11_ctd_mean,'r','linewidth',2,'displayname','CTD');
    %loglog(1/3:1/3:160,noise,'m','linewidth',2)
    %loglog(epsi_k,10.^(mmp_noise).*dTdV(1).^2,'m--','linewidth',2)
    p3 = loglog(epsi_f,A2,'Color',.8* [1 1 1],'linewidth',2,'displayname','raw(Volt^2/Hz)');
    p4 = loglog(epsi_f,AA2,'Color',.6* [1 1 1],'linewidth',2,'displayname','corrected (Volt^2/Hz)');
    %loglog(epsi_k,B,'Color',.4* [1 1 1],'linewidth',2)
    if ~isempty(indt2)
        %    loglog(epsi_k,P11_T(indt2,:),'Color',.1* [1 1 1],'linewidth',2)
    end
    
    set(gca,'XScale','log','YScale','log')
    xlabel('Hz','fontsize',20)
    ylabel('C^2/Hz','fontsize',20)
    % Title for figure
    titleStr = strrep([Meta_Data.cruise_name ' ' Meta_Data.vehicle_name ' ' Meta_Data.deployment],'_','\_');
    if isfield(Profile,'profNum')
        titleStr={sprintf('%s cast %i - temperature_ ',titleStr,Profile.profNum);...
            sprintf('dTdV = %2.2f',dTdV(2))};
    else
        titleStr={['epsitime ' datestr(Profile.epsi.dnum(1),'HH:MM:SS') ' - ' datestr(Profile.epsi.dnum(end),'HH:MM:SS')]; ...
            sprintf('dTdV = %2.2f',dTdV(2))};
    end
    title(titleStr,'fontsize',20)
    
    legend([p1,p2,p3,p4(1)],'location','northeast')
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
   
end
end

