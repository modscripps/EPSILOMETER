function [Phi,f1,noise,ax] = ttv_plot_at_tMid(obj,tMid,tscan,nSec,makeFig,saveFig,replaceData,ax)
% [Phi,f1,noise,ax] = epsiPlot_spectra_at_tMid(obj,tMid,tscan,makeFig,saveFig,replaceData,ax)
%   ** Nicole Couto copied from mod_som_calibrate_epsi_tMid.m just to start giving
%   consistent names to epsi_class plotting functions
%
%  script to calibrate the epsilometer electronics
%  Good practice: It would great to keep SOM and front end together to keep
%  track of the "system's noise"
%
%  INPUTS:
%       obj       - epsi_class object or Profile structure, inculdes
%                   Meta_Data as one of the fields
%       tMid      - midpoint of timeseries for computing spectra (seconds)
%       tScan     - length in seconds of the segment used to do an FFT (i.e tscan = NFFT * FS)
%       nSec      - length of timeseries to plot
%
% OUTPUTS:
%       Phi        - structure containing frequency spectra for each
%                    channel
%       f1         - frequency array
%       noise      - structure containing shear and fpo7 noise data
%       ax         - figure axes
%
%
%  Created by Arnaud Le Boyer on 7/28/18.
%  Copyright © 2018 Arnaud Le Boyer. All rights reserved.
%
%   March 2021 - Nicole Couto edited to take tMid AND tscan
%   April 2021 - Nicole Couto edited to include outputs
%   June  2021 - Nicole Couto edited to have ax as output and input to use
%                with epsi_realtime_spectra
% --------------------------------------------------------------------------

% If no makeFig flag, make figure by default
if nargin<7
    ax = [];
    if nargin<6
        replaceData=0;
        if nargin<5
            makeFig = 1;
            saveFig = 0;
        end
        if nargin==5
            saveFig = 1;
        end
    end
end

Meta_Data = obj.Meta_Data;
if ~isfield(obj,'plot_properties')
    obj.plot_properties = epsiSetup_set_plot_properties;
    cols = obj.plot_properties.Colors;
else
    cols = obj.plot_properties.Colors;
end
if ~any([isfield(obj,'Meta_Data'),isclassfield(obj,'Meta_Data')])
    obj.Meta_Data.mission = '';
    obj.Meta_Data.vehicle_name = '';
    obj.Meta_Data.deployment = '';
end

Meta_Data = obj.Meta_Data;
EPSI = obj.epsi;
TTV  = obj.ttv; 


% If there's no CTD data...
if ~(isfield(obj,'ttv') || isclassfield(obj,'ttv'))
    if isempty(obj.ctd)
        TTV.time_s = [];
        TTV.dtof = [];
    end
elseif isfield(obj,'ttv') || isclassfield(obj,'ttv')
    TTV = obj.ttv;
    ind0=(TTV.time_s==0);
    notNan = ~isnan(TTV.time_s);
    if sum(~ind0 & notNan)>3
        dtof_interp = interp1(TTV.time_s(~ind0 & notNan),TTV.dtof(~ind0 & notNan),EPSI.time_s);
    else
        TTV.dtof = [];
    end
end


% Get rid of nans in the data. If you're plotting in realtime, there might
% be a lot of nans at the end
EPSI = structfun(@(x) x(~isnan(x)),EPSI,'un',0);
TTV = structfun(@(x) x(~isnan(x)),TTV,'un',0);

time_s = EPSI.time_s;
time_dnum = EPSI.dnum;
%time_s = EPSI.time_s - nanmin(EPSI.time_s);
L=length(time_s);
FS=Meta_Data.AFE.FS;

% %%
% %figure('units','inch','position',[0,0,35,15]);

% plot(time_s,EPSI.t1_volt,'displayname',strrep('t1_volt','_','\_'))
% hold on
% plot(time_s,EPSI.t2_volt,'displayname','t2\_volt')
% xlabel('seconds')
% ylabel('volts')
% legend('autoupdate','off')
%
% %%
% fig2 = figure;
% plot(time_s,EPSI.s1_volt,'displayname','s1\_volt')
% hold on
% plot(time_s,EPSI.s2_volt,'displayname','s2\_volt')
% xlabel('seconds')
% ylabel('volts')
% legend
%
% %%
% fig3 = figure;
% plot(time_s,EPSI.a1_g,'displayname','a1\_g')
% hold on
% plot(time_s,EPSI.a2_g,'displayname','a2\_g')
% plot(time_s,EPSI.a3_g,'displayname','a3\_g')
% xlabel('seconds')
% ylabel('g')
% legend

% -------------------------------------------------------------------------
% Arnaud's version:
% Lseg is number of sample needed to plot 7*tscan
% if this is too long for the time serie Lseg is half the time serie (i.e. start the diag at the beginning)
% Lseg=min(FS*7*tscan,floor(L/2)-1);
% idxScan = floor(L/2)-Lseg:floor(L/2)+Lseg;

% Nicole's version:
% tscan is the number of seconds, tMid is the midpoint in time of the scan
% you want to plot (you might get this by clicking on a timeseries of
% pressure, dPdt
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
    
   

% %{ % Make y-limits based on min and max
%      maxa=max([max(a1) max(a2) max(a3)]);
%     mina=min([min(a1) min(a2) min(a3)]);
%     maxs=max([max(s1) max(s2)]);
%     mins=min([min(s1) min(s2)]);
%     maxt=max([max(t1) max(t2)]);
%     mint=min([min(t1) min(t2)]); 
% %}

% Make y-limits based on median and standard deviation
var_list = {'a1','a2','a3','s1','s2','t1','t2'};
for v=1:length(var_list)
    eval(sprintf('%smin = nanmedian(%s)-3*std(%s);',var_list{v},var_list{v},var_list{v}))
    eval(sprintf('%smax = nanmedian(%s)+3*std(%s);',var_list{v},var_list{v},var_list{v}))
end
maxa13 = max([a1max,a3max]);
mina13 = min([a1min,a3min]);
maxa2 = a2max;
mina2 = a2min;
maxs=max([s1max, s2max]);
mins=min([s1min, s2min]);
maxt=max([t1max, t2max]);
mint=min([t1min, t2min]);

    
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
    
end %end if idxScan is within data limits

%% Make the figure
% -----------------------------------------
% -----------------------------------------------
if makeFig
    
    % Set up axes
    if ~replaceData
        clear ax
        fig4 = figure;
        % Set figure size based on screen size
        defaultFigWidth = 954;
        defaultFigHeight = 954;
        screenSize = get(0,'screensize');
        mult = round(min([screenSize(3)/defaultFigWidth,screenSize(4)/defaultFigHeight]),2);
        set(fig4,'Units','pixels','Position',[1 1 defaultFigWidth*mult defaultFigHeight*mult]);
        ax(1)=subplot('Position',[0.0900    0.8969    0.8200    0.0531]);
        ax(2)=subplot('Position',[0.0900    0.8317    0.8200    0.0531]);
        ax(3)=subplot('Position',[0.0900    0.7666    0.8200    0.0531]);
        ax(4)=subplot('Position',[0.0900    0.7014    0.8200    0.0531]);
        ax(5)=subplot('Position',[0.0900    0.6363    0.8200    0.0531]);
        ax(6)=subplot('Position',[.0900    0.0500    0.8200    0.5091]);
        
        
    elseif replaceData
        if strcmp(ax(1).Tag,'tMid_spectra_a1')
            for iAx=1:length(ax)
                ax(iAx).NextPlot = 'replace';
            end
        else
            clear ax
            fig4 = figure;
            % Set figure size based on screen size
            defaultFigWidth = 954;
            defaultFigHeight = 954;
            screenSize = get(0,'screensize');
            mult = round(min([screenSize(3)/defaultFigWidth,screenSize(4)/defaultFigHeight]),2);
            set(fig4,'Units','pixels','Position',[1 1 defaultFigWidth*mult defaultFigHeight*mult]);
            ax(1)=subplot('Position',[0.0900    0.8969    0.8200    0.0531]);
            ax(2)=subplot('Position',[0.0900    0.8317    0.8200    0.0531]);
            ax(3)=subplot('Position',[0.0900    0.7666    0.8200    0.0531]);
            ax(4)=subplot('Position',[0.0900    0.7014    0.8200    0.0531]);
            ax(5)=subplot('Position',[0.0900    0.6363    0.8200    0.0531]);
            ax(6)=subplot('Position',[.0900    0.0500    0.8200    0.5091]);
            
        end
    end
    
    
    % Plot timeseries
    % --------------------
    if ~isempty(idxSeg)
        plot(ax(1),time_dnum(idxSeg),a1,'Color',cols.a1)
        hold(ax(1),'on')
        plot(ax(1),time_dnum(idxSeg),a2,'Color',cols.a2)
        %plot(ax(1),time_dnum(idxScan),squeeze(data(5,1,:)),'Color',cols.a1,'linewidth',4)
        %plot(ax(1),time_dnum(idxScan),squeeze(data(6,1,:)),'Color',cols.a2,'linewidth',4)
        hold(ax(1),'off')
        
        
        plot(ax(2),time_dnum(idxSeg),a3,'Color',cols.a3)
        
        plot(ax(3),time_dnum(idxSeg),t1,'Color',cols.t1)
        hold(ax(3),'on')
        plot(ax(3),time_dnum(idxSeg),t2,'Color',cols.t2)
        hold(ax(3),'off')
        
        plot(ax(4),time_dnum(idxSeg),s1,'Color',cols.s1)
        hold(ax(4),'on')
        plot(ax(4),time_dnum(idxSeg),s2,'Color',cols.s2)
        hold(ax(4),'off')
        
        if ~isempty(TTV.dtof)
            plot(ax(5),time_dnum(idxSeg),dtof_interp(idxSeg)*1e12,'Color',cols.T)
        end
        
        
        ax(1).YLim=[mina13 maxa13];
        ax(2).YLim=[mina2 maxa2];
        ax(3).YLim=[mint maxt];
        ax(4).YLim=[mins maxs];
        %ax(5).YLim=[-0.1 1.5];
        ylabel(ax(1),'V','FontSize',14)
        ylabel(ax(2),'V','FontSize',14)
        ylabel(ax(3),'V','FontSize',14)
        ylabel(ax(4),'V','FontSize',14)
        ylabel(ax(5),'ps','FontSize',14)
        for a=1:4
            ax(a).XTickLabel='';
            ax(a).FontSize=14;
            ax(a).XLim=[time_dnum(idxSeg(1)) time_dnum(idxSeg(end))];
        end
        ax(5).FontSize=14;
        %xlabel(ax(5),['epsitime (seconds)'],'fontsize',14)
        
        % Add datetick
        for a=1:5
            datetick(ax(a),'x','HH:MM:SS','keeplimits')
        end
        [ax(1:4).XTickLabel] = deal('');
        % shade the analyzed block
        % --------------------
        for a=1:4
            hold(ax(a),'on')
            h(a)=fill(ax(a),time_dnum(indexes), ...
                [-2 2 2 -2],[.7 .7 .7]);
            h(a).FaceAlpha=.7;
            h(a).EdgeColor = 'none';
            
            plot(ax(a),[time_dnum(idxMid) time_dnum(idxMid)], ...
                [-2 2],'k');
            plot(ax(a),[time_dnum(idxScan(1)) time_dnum(idxScan(1))], ...
                [-2 2],'k--');
            plot(ax(a),[time_dnum(idxScan(end)) time_dnum(idxScan(end))], ...
                [-2 2],'k--');
        end
        
        
        legend(ax(1),{sprintf('a1 %1.2fV',ma1),sprintf('a2 %1.2fV',ma2)},'location','northwest')
        legend(ax(2),{sprintf('a3 %1.2fV',ma3)},'location','northwest')
        legend(ax(3),{sprintf('t1 %1.2eV',mt1),sprintf('t2 %1.2eV',mt2)},'location','northwest')
        legend(ax(4),{sprintf('s1 %1.2eV',ms1),sprintf('s2 %1.2eV',ms2)},'location','northwest')
        % legend(ax(5),{'diff ramp'})
        linkaxes(ax(1:5),'x');
        
        % plot spectra
        % --------------------
        % hold(ax(6),'on')
        % l0=loglog(ax(6),f1,squeeze(nanmean(P11(1,:,:),2)),'--','Color',cmap(4,:));
        % loglog(ax(6),f1,squeeze(nanmean(P11(2,:,:),2)),'--','Color',cmap(5,:))
        % loglog(ax(6),f1,squeeze(nanmean(P11(3,:,:),2)),'--','Color',cmap(6,:))
        % loglog(ax(6),f1,squeeze(nanmean(P11(4,:,:),2)),'--','Color',cmap(7,:))
        % loglog(ax(6),f1,squeeze(nanmean(P11(5,:,:),2)),'--','Color',cmap(1,:))
        % loglog(ax(6),f1,squeeze(nanmean(P11(6,:,:),2)),'--','Color',cmap(2,:))
        % loglog(ax(6),f1,squeeze(nanmean(P11(7,:,:),2)),'--','Color',cmap(3,:))
        
        l0=loglog(ax(6),f1,squeeze(nanmean(P11bis(1,:,:),2)),'--','Color',cols.t1);
        hold(ax(6),'on')
        loglog(ax(6),f1,squeeze(nanmean(P11bis(2,:,:),2)),'--','Color',cols.t2)
        loglog(ax(6),f1,squeeze(nanmean(P11bis(3,:,:),2)),'--','Color',cols.s1)
        loglog(ax(6),f1,squeeze(nanmean(P11bis(4,:,:),2)),'--','Color',cols.s2)
        loglog(ax(6),f1,squeeze(nanmean(P11bis(5,:,:),2)),'--','Color',cols.a1)
        loglog(ax(6),f1,squeeze(nanmean(P11bis(6,:,:),2)),'--','Color',cols.a2)
        loglog(ax(6),f1,squeeze(nanmean(P11bis(7,:,:),2)),'--','Color',cols.a3)
        
        %
        l1=loglog(ax(6),f1,t1m,'d-','Color',cols.t1);
        l2=loglog(ax(6),f1,t2m,'p-','Color',cols.t2);
        l3=loglog(ax(6),f1,s1m,'d-','Color',cols.s1);
        l4=loglog(ax(6),f1,s2m,'p-','Color',cols.s2);
        l5=loglog(ax(6),f1,a1m,'d-','Color',cols.a1);
        l6=loglog(ax(6),f1,a2m,'p-','Color',cols.a2);
        l7=loglog(ax(6),f1,a3m,'s-','Color',cols.a3);
        set(ax(6),'Xscale','log','Yscale','log')
        
        % bit noise
        n20=loglog(ax(6),f1,f1*0+def_noise(20),'--','Color',[.5 .5 .5],'linewidth',2);
        n24=loglog(ax(6),f1,f1*0+def_noise(24),'--','Color',[.1 .1 .1],'linewidth',2);
        n16=loglog(ax(6),f1,f1*0+def_noise(16),'.-','Color',[.3 .3 .3],'linewidth',2);
        An1=loglog(ax(6),f1,KionixAccelnoise,'--','Color',[.1 .1 .1],'linewidth',2);
        An2=loglog(ax(6),f1,ADXLAccelnoise,'--','Color',[.1 .6 .1],'linewidth',2);
        Emp=loglog(ax(6),f1,10.^test_noise1,'m-','linewidth',2);
        Emps=loglog(ax(6),f1,10.^test_snoise1,'c-','linewidth',2);
        
        grid(ax(6),'on')
        legend([l0,l1 l2 l3 l4 l5 l6 l7 n24 n20 n16 An1 An2 Emp Emps],{'no TF','t1','t2','s1','s2','a1','a2','a3','24 bit','20 bit','16 bit','Kionix Accel noise','ADXL Accel noise','t1-noise','s1-noise'},'location','SouthWest')
        set(ax(6),'fontsize',14)
        ylabel(ax(6),'V^2 / Hz','fontsize',14)
        xlabel(ax(6),'Hz','fontsize',14)
        ax(6).XLim=[1/tscan f(end)];
        %ax(6).YLim=[0.9*def_noise(24) 10*KionixAccelnoise(1)];
        ax(6).YLim = [0.9*def_noise(24) max([t1m(:);t2m(:);s1m(:);s2m(:);a1m(:);a2m(:);a3m(:)])]; %NC changed because y-limits were never large enough
        ax(6).YLim = [0.9*def_noise(24) 1e-3];
        
        %title(ax(1),[Meta_Data.CTL.name '-' Meta_Data.CTL.rev '-' Meta_Data.CTL.SN '-' ...
            %Meta_Data.AFE.name '-' Meta_Data.AFE.rev '-' Meta_Data.AFE.SN],'fontsize',25)
        
        fig4.PaperPosition = [0 0 25 25];
        
        
        
        
        % Add tag for tracking axes if you're plotting in realtime
        ax(1).Tag = 'tMid_spectra_a1';
        
        figureStamp(getFilename)
    else
        % Add tag for tracking axes if you're plotting in realtime
        plot(ax(1),1,1);
        ax(1).Tag = 'tMid_spectra_a1';
    end %end if ~isempty(idxSeg)
    
    % Save figure
    % --------------------
    if saveFig
        img = getframe(gcf);
        imwrite(img.cdata,fullfile(Meta_Data.paths.data,['figs/epsi_' Meta_Data.deployment '_t' num2str(tMid) '.png']));
    end
    
end % end if makeFig

