function [f1,Phi,Phi_filtered,noise_curves] = epsiProcess_calc_raw_spectra(obj,tMid,tscan)
% [f1,scandata,P11,noise] = epsiProcess_calc_raw_spectra(obj,tMid,tscan)
% tMid = # seconds from the end of obj.epsi.time_s

Meta_Data = obj.Meta_Data;
EPSI = obj.epsi;

% Get rid of nans in the data. If you're plotting in realtime, there might
% be a lot of nans at the end
EPSI = structfun(@(x) x(~isnan(x)),EPSI,'un',0);

time_s = EPSI.time_s;
time_dnum = EPSI.dnum;
%time_s = EPSI.time_s - nanmin(EPSI.time_s);
L=length(time_s);
FS=Meta_Data.AFE.FS;

[~,idxMid] =  nanmin(abs(time_s - tMid));
idxScan = floor(idxMid - FS*(tscan/2)) : floor(idxMid + FS*(tscan/2));

if isempty(idxScan) || ~all(idxScan>0 & idxScan<=length(time_s))
    %If idxScan is outside data limits, do nothing
    Phi = [];
    f1 = [];
    noise = [];
    idxSeg = [];
else    

% Lscan,defined later is the length of tscan. Lseg is the length of the
    % timeseries you want to plot. Let's plot 30 seconds of data
    Lseg = FS*nSec;
    idxSeg = floor(mean(idxScan) - Lseg/2):floor(mean(idxScan) + Lseg/2);
    % If the 30-second segment goes over the length of the timeseries, or
    % begins before it, adjust accordingly
    if max(idxSeg)>L
        idxSeg = L-Lseg:L;
    elseif min(idxSeg)<1
        idxSeg = 1:Lseg;
    end
    
    % -------------------------------------------------------------------------
    
    mean_channel=@(x) (nanmean(x(idxSeg)));
    mt1=mean_channel(EPSI.t1_volt);
    mt2=mean_channel(EPSI.t2_volt);
    ms1=mean_channel(EPSI.s1_volt);
    ms2=mean_channel(EPSI.s2_volt);
    ma1=mean_channel(EPSI.a1_g);
    ma2=mean_channel(EPSI.a2_g);
    ma3=mean_channel(EPSI.a3_g);
    
    t1=detrend(EPSI.t1_volt(idxSeg));
    t2=detrend(EPSI.t2_volt(idxSeg));
    s1=detrend(EPSI.s1_volt(idxSeg));
    s2=detrend(EPSI.s2_volt(idxSeg));
    a1=detrend(EPSI.a1_g(idxSeg));
    a2=detrend(EPSI.a2_g(idxSeg));
    a3=detrend(EPSI.a3_g(idxSeg));
        

%tscan   = 5
    nb_seg=floor(time_s(end)/tscan);
    dof=5;
    if nb_seg<dof
        warning('too short. gather more data')
        %     mod_som_calibrate_epsi(Meta_Data,tscan)
    end
    
    
    % number of samples per scan (1s) in channels
    df        = 1/tscan;
    f=(df:df:FS/2)'; % frequency vector for spectra
    %% Length of the EPSI
    T       = length(EPSI.time_s);
    %% define number of scan in the EPSI
    %Lscan   = floor(tscan*FS);
    Lscan   = numel(idxScan);
    nbscan  = floor(T/Lscan);
    
    %% we compute spectra on scan with 50% overlap
    nbscan=2*nbscan-1;
    %channels=Meta_Data.PROCESS.channels;
    timeseries=Meta_Data.PROCESS.timeseries;
    nb_channels=length(timeseries);
    
    
    % NC commented: These lines of code are splitting the entire timeseries
    % into lengths of tscan and the choosing the one in the middle to plot. We
    % have already chosen the area of the timeseries to plot visually, so we
    % don't need to do all this.
    % ----
    % %% split in segment of tscan length  the time series
    % indscan = arrayfun(@(x) (1+floor(Lscan/2)*(x-1):1+floor(Lscan/2)*(x-1)+Lscan-1),1:nbscan,'un',0);
    % nb_segment=numel(indscan);
    % %% select 2 segments around the middle of time serie
    % indscan1=indscan(max(floor(nb_segment/2)-dof,1):min(floor(nb_segment/2)+dof,T));
    % idxScan=[indscan1{1}(1).*[1 1] indscan1{end}(end).*[1 1]];
    % ----
    % I'm not sure why, but giving it only one scan made get_profile_spectrum
    % crash. This seems to work.
    indscan1{1}  = idxScan;
    indscan1{2}  = idxScan;
    indscan1{3}  = idxScan;
    indexes=[indscan1{1}(1).*[1 1] indscan1{end}(end).*[1 1]];
    %idxScanLims = [idxScan(1) idxScan(end)];
    
    %% clean time series
    for c=1:nb_channels
        wh_channels=timeseries{c};
        EPSI.(wh_channels)=fillmissing(EPSI.(wh_channels)(indexes(1):indexes(3)),'linear');
        EPSI.(wh_channels)=filloutliers(EPSI.(wh_channels),'linear','movmedian',1000);
    end
    
    %% split data
    data=zeros(nb_channels,numel(indscan1),Lscan);
    for c=1:nb_channels
        wh_channels=timeseries{c};
        ind=find(cellfun(@(x) strcmp(x,wh_channels),timeseries));
        switch wh_channels
            case 't1_volt'
                data(ind,:,:) = cell2mat(cellfun(@(x) EPSI.(wh_channels)(x-indexes(1)+1),indscan1,'un',0)).';
            case 't2_volt'
                data(ind,:,:) = cell2mat(cellfun(@(x) EPSI.(wh_channels)(x-indexes(1)+1),indscan1,'un',0)).';
            case {'s1_volt','s2_volt'}
                data(ind,:,:) = cell2mat(cellfun(@(x) EPSI.(wh_channels)(x-indexes(1)+1),indscan1,'un',0)).';
            case {'a1_g','a2_g','a3_g'}
                data(ind,:,:) = cell2mat(cellfun(@(x) EPSI.(wh_channels)(x-indexes(1)+1),indscan1,'un',0)).';
        end
    end
    
    %% compute spectra
    [f1,~,P11,~]=get_profile_spectrum(data,f);
    indf1=find(f1>=0);
    indf1=indf1(1:end-1);
    f1=f1(indf1);
    P11= 2*P11(:,:,indf1);
    
    for i=1:7
        test1=detrend(squeeze(data(i,:,:)).');
        [test,f1]=pwelch(test1,[],[],4096,325);
        P11bis(i,:,:)=test.';
    end
    
    %% get TF
    h_freq=get_filters_SOM(Meta_Data,f1);
    
    %% correct the TF
    a1m=squeeze(nanmean(P11bis(5,:,:),2))./h_freq.electAccel(:);
    a2m=squeeze(nanmean(P11bis(6,:,:),2))./h_freq.electAccel(:);
    a3m=squeeze(nanmean(P11bis(7,:,:),2))./h_freq.electAccel(:);
    
    s1m=squeeze(nanmean(P11bis(3,:,:),2))./h_freq.shear(:);
    s2m=squeeze(nanmean(P11bis(4,:,:),2))./h_freq.shear(:);
    
    t1m=squeeze(nanmean(P11bis(1,:,:),2))./h_freq.electFPO7(:);
    t2m=squeeze(nanmean(P11bis(2,:,:),2))./h_freq.electFPO7(:);
    % t1m=squeeze(nanmean(P11bis(1,:,:),2));
    % t2m=squeeze(nanmean(P11bis(2,:,:),2));
    
    %% get noise
    Fn    = .5*FS;  % Nyquist frequency
    FR    = 2.5;    % Full range in Volts
    def_noise=@(x)((FR/2^x)^2 /Fn);
    KionixAccelnoise=45e-6^2+0*f1;
    ADXLAccelnoise=20e-6^2+0*f1;
    
    % polyfit the data
    logf=log10(f1);
    
    %% t1 and s1 noise
    Empnoise=log10(squeeze(nanmean(P11bis(1,:,:),2)));
    Empnoiseshear=log10(squeeze(nanmean(P11bis(3,:,:),2)));
    Emp_FPO7noise=polyfit(logf(2:end),Empnoise(2:end),3);
    Emp_shearnoise=polyfit(logf(2:end),Empnoiseshear(2:end),3);
    %test_noise=polyval(Emp_FPO7noise,logf);
    n3=Emp_FPO7noise(1);
    n2=Emp_FPO7noise(2);
    n1=Emp_FPO7noise(3);
    n0=Emp_FPO7noise(4);
    
    n3s=Emp_shearnoise(1);
    n2s=Emp_shearnoise(2);
    n1s=Emp_shearnoise(3);
    n0s=Emp_shearnoise(4);
    
    test_noise1=n0+n1.*logf+n2.*logf.^2+n3.*logf.^3;
    test_snoise1=n0s+n1s.*logf+n2s.*logf.^2+n3s.*logf.^3;
    
    noise.t1.n0 = n0;
    noise.t1.n1 = n1;
    noise.t1.n2 = n2;
    noise.t1.n3 = n3;
    
    noise.s1.n0 = n0;
    noise.s1.n1 = n1;
    noise.s1.n2 = n2;
    noise.s1.n3 = n3;
    
    %% t2 and s2 noise
    Empnoise=log10(squeeze(nanmean(P11bis(2,:,:),2)));
    Empnoiseshear=log10(squeeze(nanmean(P11bis(4,:,:),2)));
    Emp_FPO7noise=polyfit(logf(2:end),Empnoise(2:end),3);
    Emp_shearnoise=polyfit(logf(2:end),Empnoiseshear(2:end),3);
    %test_noise=polyval(Emp_FPO7noise,logf);
    n3=Emp_FPO7noise(1);
    n2=Emp_FPO7noise(2);
    n1=Emp_FPO7noise(3);
    n0=Emp_FPO7noise(4);
    
    n3s=Emp_shearnoise(1);
    n2s=Emp_shearnoise(2);
    n1s=Emp_shearnoise(3);
    n0s=Emp_shearnoise(4);
    
    test_noise2=n0+n1.*logf+n2.*logf.^2+n3.*logf.^3;
    test_snoise2=n0s+n1s.*logf+n2s.*logf.^2+n3s.*logf.^3;
    
    noise.t2.n0 = n0;
    noise.t2.n1 = n1;
    noise.t2.n2 = n2;
    noise.t2.n3 = n3;
    
    noise.s2.n0 = n0;
    noise.s2.n1 = n1;
    noise.s2.n2 = n2;
    noise.s2.n3 = n3;
    
    %% Make output structures
    Phi.t1 = squeeze(nanmean(P11bis(1,:,:),2));
    Phi.t2 = squeeze(nanmean(P11bis(2,:,:),2));
    Phi.s1 = squeeze(nanmean(P11bis(3,:,:),2));
    Phi.s2 = squeeze(nanmean(P11bis(4,:,:),2));
    Phi.a1 = squeeze(nanmean(P11bis(5,:,:),2));
    Phi.a2 = squeeze(nanmean(P11bis(6,:,:),2));
    Phi.a3 = squeeze(nanmean(P11bis(7,:,:),2));

    Phi_filtered.t1 = t1m;
    Phi_filtered.t2 = t2m;
    Phi_filtered.s1 = s1m;
    Phi_filtered.s2 = s2m;
    Phi_filtered.a1 = a1m;
    Phi_filtered.a2 = a2m;
    Phi_filtered.a3 = a3m;

    noise_curves.def_noise = def_noise;
    noise_curves.KionixAccelnoise = KionixAccelnoise;
    noise_curves.ADXLAccelnoise = ADXLAccelnoise;
    noise_curves.test_noise1 = test_noise1;
    noise_curves.test_snoise1 = test_snoise1;
end