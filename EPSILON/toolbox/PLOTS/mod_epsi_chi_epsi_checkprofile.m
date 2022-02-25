function mod_epsi_chi_epsi_checkprofile(MS,Meta_Data,l)
%function mod_epsi_chi_epsi_checkprofile(MS,Meta_Data,l)
% Feb 2019 ALB. Pass in structure MS and Meta_Data.  This will generate a
% standard movie for cast l.

close all
% l is the profile number
% l=5;
fontsize=20;
set(gcf,'Position',[500 100 2000 1500])


%
channels=Meta_Data.PROCESS.channels;

indt1=find(cellfun(@(x) strcmp(x,'t1'),channels));
indt2=find(cellfun(@(x) strcmp(x,'t2'),channels));
inds1=find(cellfun(@(x) strcmp(x,'s1'),channels));
inds2=find(cellfun(@(x) strcmp(x,'s2'),channels));
inda1=find(cellfun(@(x) strcmp(x,'a1'),channels));
inda2=find(cellfun(@(x) strcmp(x,'a2'),channels));
inda3=find(cellfun(@(x) strcmp(x,'a3'),channels));


% acceleration and coherence with shear
if ~isempty(inda1);sma1=squeeze(smoothdata(MS{l}.Pf(inda1,:,:),2,'movmean',20));
else;sma1=0.*squeeze(MS{l}.Pf(1,:,:));end
if ~isempty(inda2);sma2=squeeze(smoothdata(MS{l}.Pf(inda2,:,:),2,'movmean',20));
else;sma2=0.*squeeze(MS{l}.Pf(1,:,:));end
if ~isempty(inda3);sma3=squeeze(smoothdata(MS{l}.Pf(inda3,:,:),2,'movmean',20));
else;sma3=0.*squeeze(MS{l}.Pf(1,:,:));end
if ~isempty(inds1);sms1=squeeze(smoothdata(MS{l}.Pf(inds1,:,:),2,'movmean',20));
else;sms1=0.*squeeze(MS{l}.Pf(1,:,:));end
if ~isempty(inds2);sms2=squeeze(smoothdata(MS{l}.Pf(inds2,:,:),2,'movmean',20));
else;sms2=0.*squeeze(MS{l}.Pf(1,:,:));end



% Co11=abs(squeeze(nanmean(MS{l}.Co12(inds1,inda1-1,:,:),3)));
% Co12=abs(squeeze(nanmean(MS{l}.Co12(inds1,inda2-1,:,:),3)));
% Co13=abs(squeeze(nanmean(MS{l}.Co12(inds1,inda3-1,:,:),3)));
% Co21=abs(squeeze(nanmean(MS{l}.Co12(inds2,inda1-1,:,:),3)));
% Co22=abs(squeeze(nanmean(MS{l}.Co12(inds2,inda2-1,:,:),3)));
% Co23=abs(squeeze(nanmean(MS{l}.Co12(inds2,inda3-1,:,:),3)));
if ~isempty(inda1);Co11=abs(squeeze(smoothdata(MS{l}.Co12(inds1,inda1-1,:,:),3,'movmean',20)));
else Co11=0.*squeeze(MS{l}.Pf(1,:,:));end
if ~isempty(inda2);Co12=abs(squeeze(smoothdata(MS{l}.Co12(inds1,inda2-1,:,:),3,'movmean',20)));
else Co12=0.*squeeze(MS{l}.Pf(1,:,:));end
if ~isempty(inda3);Co13=abs(squeeze(smoothdata(MS{l}.Co12(inds1,inda3-1,:,:),3,'movmean',20)));
else Co13=0.*squeeze(MS{l}.Pf(1,:,:));end
if ~isempty(inda1);Co21=abs(squeeze(smoothdata(MS{l}.Co12(inds2,inda1-1,:,:),3,'movmean',20)));
else Co21=0.*squeeze(MS{l}.Pf(1,:,:));end
if ~isempty(inda2);Co22=abs(squeeze(smoothdata(MS{l}.Co12(inds2,inda2-1,:,:),3,'movmean',20)));
else Co22=0.*squeeze(MS{l}.Pf(1,:,:));end
if ~isempty(inda3);Co23=abs(squeeze(smoothdata(MS{l}.Co12(inds2,inda3-1,:,:),3,'movmean',20)));
else Co23=0.*squeeze(MS{l}.Pf(1,:,:));end

% if isempty(Co11);Co11=zeros(1,size(MS{l}.Pf,3));end
% if isempty(Co12);Co12=zeros(1,size(MS{l}.Pf,3));end
% if isempty(Co13);Co13=zeros(1,size(MS{l}.Pf,3));end
% if isempty(Co21);Co21=zeros(1,size(MS{l}.Pf,3));end
% if isempty(Co22);Co22=zeros(1,size(MS{l}.Pf,3));end
% if isempty(Co23);Co23=zeros(1,size(MS{l}.Pf,3));end




dTdV=[Meta_Data.epsi.t1.dTdV Meta_Data.epsi.t2.dTdV];

% noise floor
logf=log10(MS{l}.f);
h_freq=get_filters_MADRE(Meta_Data,MS{l}.f);
switch Meta_Data.MAP.temperature
    case 'Tdiff'
        FPO7noise=load(fullfile(Meta_Data.paths.calibration,'FPO7_noise.mat'),'n0','n1','n2','n3');
    otherwise
        FPO7noise=load(fullfile(Meta_Data.paths.calibration,'FPO7_notdiffnoise.mat'),'n0','n1','n2','n3');
end
n0=FPO7noise.n0; n1=FPO7noise.n1; n2=FPO7noise.n2; n3=FPO7noise.n3;
tnoise=10.^(n0+n1.*logf+n2.*logf.^2+n3.*logf.^3);


shearnoise=load(fullfile(Meta_Data.paths.calibration,'shear_noise.mat'),'n0s','n1s','n2s','n3s');
n0s=shearnoise.n0s; n1s=shearnoise.n1s; n2s=shearnoise.n2s; n3s=shearnoise.n3s;
snoise=10.^(n0s+n1s.*logf+n2s.*logf.^2+n3s.*logf.^3);


% % movie stuff
v = VideoWriter(sprintf('%s_cast%i%s',Meta_Data.deployment,l,'.avi'));
v.FrameRate=5;
open(v)

% TG spectra 
ax(1)=axes('Position',[.05 .06 .45 .6]);
% chi profiles
ax(2)=axes('Position',[.066 .63 .08 .3]);
% epsi profiles
ax(3)=axes('Position',[.15 .63 .08 .3]);
% T-S profiles
ax(4)=axes('Position',[.235 .63 .08 .3]);
% w profile
ax(5)=axes('Position',[.32 .63 .08 .3]);
%  shear spectra 
ax(6)=axes('Position',[.51 .06 .44 .6]);
%  accel spectra 
ax(7)=axes('Position',[.7 .83 .24 .15]);
%  coherence spectra 
ax(8)=axes('Position',[.7 .67 .24 .15]);

a=2;
plot(ax(a),MS{l}.chi,MS{l}.pr)
a=3;
plot(ax(a),MS{l}.epsilon,MS{l}.pr)
a=4;
[ax1,hl1,hl2]=plotxx(MS{l}.t,MS{l}.pr,MS{l}.s,MS{l}.pr,{'',''},{'',''},ax(a));
a=5;
plot(ax(a),MS{l}.w,MS{l}.pr)


annotation('textbox',...
    [.56 .665 .135 .305],...
    'String',{Meta_Data.mission,...
    Meta_Data.vehicle_name,...
    Meta_Data.deployment,...
    Meta_Data.path_mission,...
    ['MADRE ' Meta_Data.MADRE.SN ' rev ' Meta_Data.MADRE.rev],...
    ['MAP '   Meta_Data.MAP.SN   ' rev ' Meta_Data.MAP.rev],...
    [Meta_Data.aux1.name ' ' Meta_Data.aux1.SN],...
    sprintf('s1 - SN: %s - Sv: %3.2f', Meta_Data.epsi.s1.SN,Meta_Data.epsi.s1.Sv),...
    sprintf('     %s - %s', Meta_Data.epsi.s1.ADCfilter,Meta_Data.epsi.s1.ADCconf),...
    sprintf('s2 - SN: %s - Sv:%3.2f ',Meta_Data.epsi.s2.SN,Meta_Data.epsi.s2.Sv),...
    sprintf('    %s - %s', Meta_Data.epsi.s2.ADCfilter,Meta_Data.epsi.s2.ADCconf),...
    sprintf('t1 - SN: %s - dTdV: %3.2f ',Meta_Data.epsi.t1.SN,Meta_Data.epsi.t1.dTdV),...
    sprintf('    %s - %s',Meta_Data.epsi.t1.ADCfilter,Meta_Data.epsi.t1.ADCconf),...
    sprintf('t2 - SN: %s - dTdV: %3.2f',Meta_Data.epsi.t2.SN,Meta_Data.epsi.t2.dTdV),...
    sprintf('    %s - %s',Meta_Data.epsi.t2.ADCfilter,Meta_Data.epsi.t2.ADCconf),...
    },...
    'FontSize',14,...
    'FontName','Arial',...
    'LineStyle','-',...
    'EdgeColor','k',...
    'LineWidth',2,...
    'BackgroundColor',[0.9  0.9 0.9],...
    'Color','k');


dTdV=[Meta_Data.epsi.t1.dTdV,Meta_Data.epsi.t2.dTdV];
Sv=[Meta_Data.epsi.s1.Sv,Meta_Data.epsi.s2.Sv];
Gr=9.81;
for k=1:2:length(MS{l}.kvis)
    
    
    % coherence plot
    a=8;
    H=semilogx(ax(a),MS{l}.f,Co11(k,:),'b');
    hold(ax(a),'on')
    I=semilogx(ax(a),MS{l}.f,Co12(k,:),'r');
    J=semilogx(ax(a),MS{l}.f,Co13(k,:),'k');
    grid(ax(a),'on')
    legend(ax(a),'a1','a2','a3','location','northwest')
    ax(a).YAxisLocation='right';
    xlim(ax(a),MS{l}.f([1 end]))
    ylim(ax(a),[0 1])
    ax(a).FontSize=fontsize;
    ax(a).XScale='log';
    ax(a).YScale='linear';
    ylabel(ax(a),'Coherence s1','fontsize',fontsize)
    xlabel(ax(a),'Hz','fontsize',fontsize)
    set(ax(a),'Xtick',[1 10 100])
    %acceleration plot
    a=7;
    K=loglog(ax(a),MS{l}.f,sma1(k,:),'b');
    hold(ax(a),'on')
    L=loglog(ax(a),MS{l}.f,sma2(k,:),'r');
    M=loglog(ax(a),MS{l}.f,sma3(k,:),'k');
    N=loglog(ax(a),MS{l}.f,sms1(k,:),'m');
    set(ax(a),'Xticklabel','')
    grid(ax(a),'on')
    legend(ax(a),'a1','a2','a3','s1','location','southwest')
    set(ax(a),'Xticklabel','')
    ax(a).YAxisLocation='right';
    xlim(ax(a),MS{l}.f([1 end]))
    ylim(ax(a),[1e-10 1e-3])
    ax(a).XScale='log';
    ax(a).YScale='log';
    ax(a).FontSize=fontsize;
    set(ax(a),'Xtick',[1 10 100])
    ylabel(ax(a),'g^2/Hz','fontsize',fontsize)
    
    
    % noise stuff
    k_noise=MS{l}.f./MS{l}.w(k);
    noise_t=tnoise.*dTdV(1).^2./h_freq.FPO7(MS{l}.w(k));
    tnoise_k= (2*pi*k_noise).^2 .* noise_t.*MS{l}.w(k);        % T1_k spec  as function of k
    

    TFshear=(Sv(1).*MS{l}.w(k)/(2*Gr)).^2 .* h_freq.shear.* haf_oakey(MS{l}.f,MS{l}.w(k));
    snoise_k= (2*pi*k_noise).^2 .* snoise.*MS{l}.w(k)./TFshear;        % T1_k spec  as function of k

    
    [kbatch11,Pbatch11] = batchelor(MS{l}.epsilon(k,1),MS{l}.chi(k,1), ...
        MS{l}.kvis(k),MS{l}.ktemp(k));
    [kbatch12,Pbatch12] = batchelor(MS{l}.epsilon(k,1),MS{l}.chi(k,2), ...
        MS{l}.kvis(k),MS{l}.ktemp(k));
    [kbatch21,Pbatch21] = batchelor(MS{l}.epsilon(k,2),MS{l}.chi(k,1), ...
        MS{l}.kvis(k),MS{l}.ktemp(k));
    [kbatch22,Pbatch22] = batchelor(MS{l}.epsilon(k,2),MS{l}.chi(k,2), ...
        MS{l}.kvis(k),MS{l}.ktemp(k));

    smTG1=smoothdata(MS{l}.PphiT_k(k,:,1),'movmean',10);
    smTG2=smoothdata(MS{l}.PphiT_k(k,:,2),'movmean',10);
    
    loglog(ax(1),MS{l}.k,MS{l}.PphiT_k(k,:,1),'--','Color',.8*[.5 1 .5])
    hold(ax(1),'on')
    loglog(ax(1),MS{l}.k,smTG1,'Color',.8*[.5 1 .7],'linewidth',2)
    loglog(ax(1),MS{l}.k,MS{l}.PphiT_k(k,:,2),'--','Color',.7*[1 1 1])
    loglog(ax(1),MS{l}.k,smTG2,'Color',.3*[1 1 1],'linewidth',2)
    loglog(ax(1),k_noise,tnoise_k,'c','linewidth',1)
    if ~isempty(indt1)
        indkc=find(MS{l}.k>MS{l}.kcfpo7(k,1),1,'first');
    scatter(ax(1),MS{l}.kcfpo7(k,1),smTG1(indkc),500,.8*[.5 1 .7],'filled','d','MarkerEdgeColor','y','linewidth',3)
    end
    if ~isempty(indt2)
        indkc=find(MS{l}.k>MS{l}.kcfpo7(k,1),1,'first');
    scatter(ax(1),MS{l}.kcfpo7(k,2),smTG2(indkc),500,.3*[1 1 1],'filled','p','MarkerEdgeColor','y','linewidth',3)
    end
    
    lolo=loglog(ax(1),kbatch11,Pbatch11,'Color',[1 0 1]);
    loglog(ax(1),kbatch12,Pbatch12,'Color',.75*[1 0 1])
    loglog(ax(1),kbatch21,Pbatch21,'Color',.5*[1 0 1])
    loglog(ax(1),kbatch22,Pbatch22,'Color',.3*[1 0 1])
    hold(ax(1),'off')
    set(ax(1),'Xscale','log','Yscale','log')
    set(ax(1),'fontsize',fontsize)
    legend(ax(1),'t1','t1smooth','t2','t2smooth','T-noise',...
                 't1_{cutoff}','t2_{cutoff}','batch11','batch12','batch21','batch22',...
                 'location','southwest' )
    xlim(ax(1),[6e-1 400])
    ylim(ax(1),[1e-10 1e-1])
    grid(ax(1),'on')
    xlabel(ax(1),'k (cpm)','fontsize',fontsize)
    ylabel(ax(1),'\phi^2_{TG} (C^2 m^{-2} / cpm)','fontsize',fontsize)
    %title(ax(1),,'position',[30 0.2]) % postion in x-y units 

    %plot chi
    a=2;
%    plot(ax(a),MS{l}.chi,MS{l}.pr)
    hold(ax(a),'on')
    ax(a).YDir='reverse';
    A=scatter(ax(a),MS{l}.chi(k,1),MS{l}.pr(k),100,'k','d','filled');
    B=scatter(ax(a),MS{l}.chi(k,2),MS{l}.pr(k),100,'k','p','filled');
    hold(ax(a),'off')
    legend(ax(a),'t1','t2','location','northeast')
    set(ax(a),'Xscale','log','Yscale','linear')
    set(ax(a),'Xtick',[1e-12 1e-10 1e-8 1e-6 1e-4])
    xlim(ax(a),[1e-12 max(MS{l}.chi(:))])
    ylim(ax(a),[min(MS{l}.pr) max(MS{l}.pr)])
    set(ax(a),'fontsize',15)
    ylabel(ax(a),'Depth (m)','fontsize',fontsize)
    xlabel(ax(a),'\chi (K^2 s^{-1}) ','fontsize',fontsize)

    %plot epsilon
    a=3;
    %plot(ax(a),MS{l}.epsilon,MS{l}.pr)
    hold(ax(a),'on')
    ax(a).YDir='reverse';
    C=scatter(ax(a),MS{l}.epsilon(k,1),MS{l}.pr(k),100,'k','d','filled');
    D=scatter(ax(a),MS{l}.epsilon(k,2),MS{l}.pr(k),100,'k','p','filled');
    hold(ax(a),'off')
    legend(ax(a),'s1','s2','location','northeast')
    set(ax(a),'Xscale','log','Yscale','linear')
    set(ax(a),'Xtick',[ 1e-10 1e-8 1e-6 1e-4])
    xlim(ax(a),[1e-12 max(MS{l}.epsilon(:))])
    ylim(ax(a),[min(MS{l}.pr) max(MS{l}.pr)])
    set(ax(a),'fontsize',15)
    set(ax(a),'Yticklabel','')
    xlabel(ax(a),'\epsilon (W kg^{-1}) ','fontsize',fontsize)

    %plot temperature
    a=4;
    %plot(ax(a),MS{l}.t,MS{l}.pr)
    %axes(ax(a));
    %[ax1,hl1,hl2]=plotxx(MS{l}.t,MS{l}.pr,MS{l}.s,MS{l}.pr,{'',''},{'',''});
    %hold(ax(1),'on')
    Dt=(max(MS{l}.t)-min(MS{l}.t))./4; 
    hold(ax1(1),'on')
    ax1(1).YDir='reverse';
    E=scatter(ax1(1),MS{l}.t(k),MS{l}.pr(k),100,'k','d','filled');
    hold(ax1(1),'off')
    hold(ax1(2),'on')
    F=scatter(ax1(2),MS{l}.s(k),MS{l}.pr(k),100,'k','p','filled');
    %hold(ax(1),'off')
    hold(ax1(2),'off')
    %legend(ax(a),'T','location','northeast')
    set(ax1(1),'Xscale','linear','Yscale','linear')
    set(ax1(1),'Xtick',floor(min(MS{l}.t)):2:floor(max(MS{l}.t)))
    xlim(ax1(1),[min(MS{l}.t) max(MS{l}.t)])
    ylim(ax1(1),[min(MS{l}.pr) max(MS{l}.pr)])
    xlim(ax1(2),[min(MS{l}.s) max(MS{l}.s)])
    ylim(ax1(2),[min(MS{l}.pr) max(MS{l}.pr)])
    set(ax1(1),'Xtick',min(MS{l}.t)+Dt:Dt:max(MS{l}.t)-Dt)
    set(ax1(1),'XtickLabel',num2str((min(MS{l}.t)+Dt:Dt:max(MS{l}.t)-Dt).','%2.1f'))
    ax1(1).XTickLabelRotation=25;
    set(ax1(1),'fontsize',15)
    set(ax1(1),'Yticklabel','')
    xlabel(ax1(1),'SBE T (C) ','fontsize',fontsize)
    
    Ds=(max(MS{l}.s)-min(MS{l}.s))./4;
    set(ax1(2),'Xscale','linear','Yscale','linear')
    set(ax1(2),'Xtick',min(MS{l}.s)+Ds:Ds:max(MS{l}.s)-Ds)
    set(ax1(2),'XtickLabel',num2str((min(MS{l}.s)+Ds:Ds:max(MS{l}.s)-Ds).','%2.1f'))
    ax1(2).XTickLabelRotation=25;
    xlim(ax1(2),[min(MS{l}.s) max(MS{l}.s)])
    ylim(ax1(2),[min(MS{l}.pr) max(MS{l}.pr)])
    set(ax1(2),'fontsize',15)
    set(ax1(2),'Yticklabel','')
    ax1(2).YDir='reverse';
    xlabel(ax1(2),'SBE S (psu) ','fontsize',fontsize)

    %plot speed
    a=5;
    %plot(ax(a),MS{l}.w,MS{l}.pr)
    hold(ax(a),'on')
    ax(a).YDir='reverse';
    G=scatter(ax(a),MS{l}.w(k),MS{l}.pr(k),100,'k','d','filled');
    hold(ax(a),'off')
    legend(ax(a),'w','location','northeast')
    set(ax(a),'Xscale','linear','Yscale','linear')
    set(ax(a),'Xtick',[.4  .6 .8 1])
    xlim(ax(a),[ .5*nanmean(MS{l}.w) 1.5*nanmean(MS{l}.w)])
    ylim(ax(a),[min(MS{l}.pr) max(MS{l}.pr)])
    set(ax(a),'fontsize',15)
    set(ax(a),'Yticklabel','')
    xlabel(ax(a),'speed (m s^{-1}) ','fontsize',fontsize)

    
    %plot shear
    a=6;
    [kpan1,Ppan1] = panchev(MS{l}.epsilon(k,1),MS{l}.kvis(k));
    [kpan2,Ppan2] = panchev(MS{l}.epsilon(k,2),MS{l}.kvis(k));
    smS1=smoothdata(MS{l}.Pshear_k(k,:,1),'movmean',10);
    smS2=smoothdata(MS{l}.Pshear_k(k,:,2),'movmean',10);
    kcindex1=find(MS{l}.k<MS{l}.kc(k,1),1,'last');
    kcindex2=find(MS{l}.k<MS{l}.kc(k,2),1,'last');
    
    loglog(ax(a),MS{l}.k,MS{l}.Pshear_k(k,:,1),'--','Color',.8*[.5 1 .5])
    hold(ax(a),'on')
    loglog(ax(a),MS{l}.k,smS1,'Color',.8*[.5 1 .7],'linewidth',2)
    loglog(ax(a),MS{l}.k,MS{l}.Pshear_k(k,:,2),'--','Color',.7*[1 1 1])
    loglog(ax(a),MS{l}.k,smS2,'Color',.3*[1 1 1],'linewidth',2)
    loglog(ax(a),k_noise,snoise_k,'c','linewidth',1)
    scatter(ax(a),MS{l}.k(kcindex1),smS1(kcindex1),500,.8*[.5 1 .7],'filled','d','MarkerEdgeColor','c','linewidth',2)
    scatter(ax(a),MS{l}.k(kcindex2),smS2(kcindex2),500,.3*[1 1 1],'filled','p','MarkerEdgeColor','c','linewidth',2)
    
    loglog(ax(a),kpan1,Ppan1,'Color',.75*[1 0 1])
    loglog(ax(a),kpan2,Ppan2,'Color',.4*[1 0 1])
    hold(ax(a),'off')
    set(ax(a),'Xscale','log','Yscale','log')
    set(ax(a),'fontsize',fontsize)
    legend(ax(a),'s1','s1smooth','s2','s2smooth','noise','s1_{cutoff}','s2_{cutoff}','Panchev1','Panchev2','location','southwest')
    xlim(ax(a),[6e-1 400])
    ylim(ax(a),[1e-10 1e-1])
    grid(ax(a),'on')
    xlabel(ax(a),'k (cpm)','fontsize',fontsize)
    ylabel(ax(a),'\phi^2_{shear} (s^{-2} / cpm)','fontsize',fontsize)
    ax(6).YAxisLocation='right';

    
    
    annotation('textbox',...
    [.41 .7 .14 .27],...
    'String',{datestr(MS{l}.time(k)),...
              sprintf('pressure=%3.1f db',MS{l}.pr(k)),...
              sprintf('speed=%1.2f m/s',MS{l}.w(k)),...
              sprintf('temperature=%2.2f (C)',MS{l}.t(k)),...
              sprintf('salinity=%2.2f (psu)',MS{l}.s(k)),...
              sprintf('kinematic viscosity =%1.1e m^2 s^{-1}',MS{l}.kvis(k)),...
              sprintf('scalar diffusivity =%1.1e m^2 s^{-1}',MS{l}.ktemp(k)),' ',...
              ['\epsilon_{1,2}' sprintf('=%1.2e, %1.2e (W kg^{-1})',MS{l}.epsilon(k,1),MS{l}.epsilon(k,2))],...
              ['\chi_{1,2}'     sprintf('=%1.2e, %1.2e (K^2 s^{-1})',MS{l}.chi(k,1),MS{l}.chi(k,2))],...
              },...
    'FontSize',14,...
    'FontName','Arial',...
    'LineStyle','-',...
    'EdgeColor','k',...
    'LineWidth',2,...
    'BackgroundColor',[0.9  0.9 0.9],...
    'Color','k');


    pause(.001)
    % movie stuff
    frame=getframe(gcf);
    writeVideo(v,frame)
    if k<length(MS{l}.kvis)
        delete(A);
        delete(B);
        delete(C);
        delete(D);
        delete(E);
        delete(F);
        delete(G);
        delete(H);
        delete(I);
        delete(J);
        delete(K);
        delete(L);
        delete(M);
        delete(N);
    end
   
end
% movie stuff
close(v)
