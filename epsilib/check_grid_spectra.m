function check_grid_spectra()



gridfilename='griddedProfiles.mat';
gridfilepath= fullfile('profiles',gridfilename);
if isfile(gridfilepath)
    load(gridfilepath,'GRID');
else
    warning(['Can not find ' gridfilepath])
    return
end


%% clean up data
nanmask=~isnan(GRID.dnum);
sgth=real(GRID.sgth);
sgth=filloutliers(sgth.','linear','movmean',100).';


%% plot GRID
close all
figgrid=figure('units','centimeters','Position',[0 40 20 40]);

axG(1)=subplot('Position',[.1 0.79 0.8 .15]); %epsilon1
axG(2)=subplot('Position',[.1 0.56 0.8 .15]); %epsilon2
axG(3)=subplot('Position',[.1 0.39 0.8 .15]); %chi1
axG(4)=subplot('Position',[.1 0.22 0.8 .15]); %chi2
axG(5)=subplot('Position',[.1 0.05 0.8 .15]); %chi2


a=1;
pcolor(axG(a),GRID.dnum,GRID.pr,log10(real(GRID.epsilon_final)));
caxis(axG(a),[-10 -8])
hold(axG(a),'on')
contour(axG(a),GRID.dnum(nanmask),GRID.pr,sgth(:,nanmask),20,'k')
hold(axG(a),'off')
title(axG(a),'epsilon final','FontSize',20,'FontName','time new roman')

a=2;
pcolor(axG(a),GRID.dnum,GRID.pr,log10(real(GRID.epsilon_co1)));
caxis(axG(a),[-10 -8])
hold(axG(a),'on')
contour(axG(a),GRID.dnum(nanmask),GRID.pr,sgth(:,nanmask),20,'k')
hold(axG(a),'off')
title(axG(a),'epsilon co 1','FontSize',20,'FontName','time new roman')


a=3;
pcolor(axG(a),GRID.dnum,GRID.pr,log10(real(GRID.epsilon_co2)));
caxis(axG(a),[-10 -8])
hold(axG(a),'on')
contour(axG(a),GRID.dnum(nanmask),GRID.pr,sgth(:,nanmask),20,'k')
hold(axG(a),'off')
title(axG(a),'epsilon co 2','FontSize',20,'FontName','time new roman')


a=4;
pcolor(axG(a),GRID.dnum,GRID.pr,log10(real(GRID.chi1)));
caxis(axG(a),[-10 -6])
hold(axG(a),'on')
contour(axG(a),GRID.dnum(nanmask),GRID.pr,sgth(:,nanmask),20,'k')
hold(axG(a),'off')
title(axG(a),'chi 1','FontSize',20,'FontName','time new roman')

a=5;
pcolor(axG(a),GRID.dnum,GRID.pr,log10(real(GRID.chi2)));
caxis(axG(a),[-10 -6])
hold(axG(a),'on')
contour(axG(a),GRID.dnum(nanmask),GRID.pr,sgth(:,nanmask),20,'k')
hold(axG(a),'off')
title(axG(a),'chi 2','FontSize',20,'FontName','time new roman')

colorbar(axG(3),'Position',[.92 0.54 0.01 .4])
colorbar(axG(5),'Position',[.92 0.08  0.01 .3])

for a=1:length(axG)
    axG(a).FontName='time new roman';
    axG(a).FontSize=15;
    shading(axG(a),'flat');
    axG(a).XTick=GRID.dnum(1):1/24:GRID.dnum(end);    
    axis(axG(a),'ij')
    ylabel(axG(a),'Depth','FontSize',20,'FontName','time new roman')
    if a==length(axG)
        timeaxis=GRID.dnum(~isnan(GRID.dnum));
        timeaxis=timeaxis(1):1/24:timeaxis(end);
        axG(a).XTickLabel=datestr(timeaxis.','HH-MM');
        axG(a).XTickLabelRotation=45;
        xlabel(axG(a),datestr(timeaxis,'dd-mm-yyyy'),'FontSize',20,'FontName','time new roman')
    else
        axG(a).XTickLabels='';
    end
end


%%

plot_flag=1;
first_time=1;
while plot_flag
    disp("Click on the epsilon profile.")
    if first_time==1
        idx=ginput(1);
        first_time=0;
    end

    select_dnum=abs(GRID.dnum-idx(1));
    id_profile=find(select_dnum==min(select_dnum));
    select_depth=idx(2);

    filename=sprintf('Profile%03i.mat',id_profile);
    rootpath=strsplit(pwd,'/');
    filepath= ...
        fullfile('profiles',filename);
    if isfile(filepath)
        load(filepath,'Profile');
    else
        warning(['Can not find ' filepath])
        return
    end


    idx1=find(Profile.pr>select_depth,1,'first');
    idx_ctd=find(Profile.ctd.P>select_depth,1,'first');



    %highlight profile and position on the GRID plots
    for a=1:4
        hold(axG(a),'on')
        la1(a)=plot(axG(a),GRID.dnum(id_profile).* [1 1],[GRID.z([1 end])],'m','LineWidth',2);
        sca1(a)=scatter(axG(a),GRID.dnum(id_profile),select_depth,100,'dy','filled');
        hold(axG(a),'on')
    end

    %% plot profile and spectra


    figure('units','centimeters','Position',[25 40 60 40])
    %
    strfontsize=15;
    % function plot_epsi_spectra(Profile)
    try
        T=filloutliers(real(Profile.ctd.T),'linear','movmean',floor(length(Profile.ctd.T)/10));
        S=filloutliers(real(Profile.ctd.S),'linear','movmean',floor(length(Profile.ctd.T)/10));
    catch
        T=real(Profile.ctd.T);
        S=real(Profile.ctd.S);
    end
    % T=Profile.ctd.T;
    % S=Profile.ctd.S;

    if isfield(Profile,'vnav')
        Profile.vnav.pr=interp1(Profile.ctd.dnum,Profile.ctd.P,Profile.vnav.dnum);
        idx_vnav=find(Profile.vnav.pr>select_depth,1,'first');
        Profile.vnav=vnav_pitch_roll_heading(Profile.vnav);
    else
        Profile.vnav.pitch   = Profile.pr.* nan;
        Profile.vnav.roll    = Profile.pr.* nan;
        Profile.vnav.heading = Profile.pr.* nan;
        Profile.vnav.dnum    = Profile.dnum;
        Profile.vnav.pr      = Profile.pr;
    end

    a=1
    ax(a)=subplot('Position',[.04 .08 .08 .8]);
    semilogx(ax(a),Profile.epsilon_co,Profile.pr);
    hold(ax(a),"on")
    semilogx(ax(a),Profile.epsilon_final,Profile.pr,'k','linewidth',2);
    hold(ax(a),"off")
    leg(a)=legend(ax(a),'s1','s2','location','northoutside');
    leg(a).Position=[.07 .95 .03 .02];
    ax(a).FontSize=strfontsize;
    ax(a).FontName='Time New Roman';
    grid(ax(a),'on')
    axis(ax(a),'ij')
    ax(a).XLim=[min(real(Profile.epsilon_co(:))) max(real(Profile.epsilon_co(:)))];
    ax(a).YLim=[max(min(Profile.pr),0) max(Profile.pr)];
    xlabel(ax(a),'\epsilon','fontsize',strfontsize,'fontname','time new roman')
    ylabel(ax(a),'Depth [m]','fontsize',strfontsize,'fontname','time new roman')

    a=a+1;
    ax(a)=subplot('Position',[.13 .08 .08 .8]);
    semilogx(ax(a),Profile.chi,Profile.pr);
    leg(a)=legend(ax(a),'t1','t2','location','northoutside');
    leg(a).Position=[.17 .95 .03 .02];
    ax(a).FontSize=strfontsize;
    ax(a).FontName='Time New Roman';
    grid(ax(a),'on')
    axis(ax(a),'ij')
    ax(a).XLim=[min(Profile.chi(:)) max(Profile.chi(:))];
    ax(a).YLim=[max(min(Profile.pr),0) max(Profile.pr)];
    ax(a).YTickLabel='';
    xlabel(ax(a),'\chi','fontsize',strfontsize,'fontname','time new roman')

    a=a+1;
    ax(a)=subplot('Position',[.22 .08 .05 .8]);
    [ax1,h(1),h(2)]=plotxx(T,Profile.ctd.P, ...
        S,Profile.ctd.P, ...
        {'T [˚C]','S [psu]'},{'',''},ax(a));
    axis(ax(a),'ij')
    ax(a).FontSize=strfontsize;
    ax(a).FontName='Time New Roman';
    grid(ax(a),'on')
    axis(ax(a),'ij')
    axis(ax1,'ij')
    ax(a).XLim=[min(T) max(T)];
    ax1(2).XLim=[min(S) max(S)];
    ax(a).YLim=[max(min(Profile.pr),0) max(Profile.pr)];
    ax1(2).YLim=[max(min(Profile.pr),0) max(Profile.pr)];
    ax(a).YTickLabel='';
    ax1(2).YTickLabel='';
    xlabel(ax(a),'T','fontsize',strfontsize,'fontname','time new roman')
    xlabel(ax1(2),'S','fontsize',strfontsize,'fontname','time new roman')
    ax1(2).FontSize=strfontsize;


    a=a+1;
    ax(a)=subplot('Position',[.4 .08 .05 .8]);
    plot(ax(a),Profile.w,Profile.pr);
    axis(ax(a),'ij')
    ax(a).FontSize=strfontsize;
    ax(a).FontName='Time New Roman';
    grid(ax(a),'on')
    axis(ax(a),'ij')
    ax(a).XLim=[min(Profile.w) max(Profile.w)];
    ax(a).YLim=[max(min(Profile.pr),0) max(Profile.pr)];
    ax(a).YTickLabel='';
    xlabel(ax(a),'w [m/s]','fontsize',strfontsize,'fontname','time new roman')

    a=a+1;
    ax(a)=subplot('Position',[.28 .08 .05 .8]);
    if isfield(Profile,'vnav')
        plot(ax(a),Profile.vnav.pitch,Profile.vnav.pr);
        axis(ax(a),'ij')
        ax(a).FontSize=strfontsize;
        ax(a).FontName='Time New Roman';
        grid(ax(a),'on')
        axis(ax(a),'ij')
        ax(a).XLim=[min(Profile.vnav.pitch) max(Profile.vnav.pitch)];
        ax(a).YLim=[max(min(Profile.pr),0) max(Profile.pr)];
        ax(a).YTickLabel='';
        xlabel(ax(a),'pitch','fontsize',strfontsize,'fontname','time new roman')
    end

    a=a+1;
    ax(a)=subplot('Position',[.34 .08 .05 .8]);
    axis(ax(a),'ij')
    ax(a).FontSize=strfontsize;
    ax(a).FontName='Time New Roman';
    grid(ax(a),'on')
    str_title=sprintf('%s %s Profile %i ',            ...
        Profile.Meta_Data.mission,      ...
        Profile.Meta_Data.vehicle_name, ...
        Profile.profNum);
    if isfield('vnav',Profile)
        plot(ax(a),unwrap(Profile.vnav.roll)+pi,Profile.vnav.pr);
        text(ax(a),-pi/2,Profile.pr(1)-15,{rootpath{end},...
            str_title},...
            'fontsize',strfontsize,...
            'fontname','time new roman',...
            'Interpreter','none')

        grid(ax(a),'on')
        axis(ax(a),'ij')
        ax(a).XLim=[min(unwrap(Profile.vnav.roll)+pi) max(unwrap(Profile.vnav.roll)+pi)];
        ax(a).YLim=[max(min(Profile.pr),0) max(Profile.pr)];
        ax(a).YTickLabel='';
        xlabel(ax(a),'roll','fontsize',strfontsize,'fontname','time new roman')
    else
        plot(ax(a),Profile.pr.*0,Profile.pr);
        text(ax(a),-2,Profile.pr(1)-3,{rootpath{end},...
            str_title},...
            'fontsize',strfontsize,...
            'fontname','time new roman',...
            'Interpreter','none')

        axis(ax(a),'ij')
        ax(a).XLim=[-1 1];
        ax(a).YLim=[max(min(Profile.pr),0) max(Profile.pr)];
    end




    a=a+1;
    ax(a)=subplot('Position',[.53 .5 .45 .13]);
    title(ax(a),'Shear','fontsize',strfontsize,'fontname','time new roman')
    ylabel(ax(a),'s^{-2}/cpm','fontsize',strfontsize,'fontname','time new roman')
    grid(ax(a),'on')
    ax(a).FontSize=20;
    ax(a).XTickLabel='';
    ax(a).FontName='Time New Roman';

    a=a+1;
    ax(a)=subplot('Position',[.53 .35 .45 .13]);
    title(ax(a),'FPO7')
    ylabel(ax(a),'˚K^2 m^{-2}/cpm','fontsize',20,'fontname','time new roman')
    grid(ax(a),'on')
    ax(a).XTickLabel='';
    ax(a).FontSize=strfontsize;
    ax(a).FontName='Time New Roman';

    a=a+1;
    ax(a)=subplot('Position',[.53 .21 .45 .13]);
    title(ax(a),'Acceleration')
    ylabel(ax(a),' g^2 /cpm','fontsize',strfontsize,'fontname','time new roman')
    ax(a).XTickLabel='';
    grid(ax(a),'on')
    ax(a).FontSize=strfontsize;
    ax(a).FontName='Time New Roman';

    a=a+1;
    ax(a)=subplot('Position',[.53 .08 .45 .08]);
    grid(ax(a),'on')
    title(ax(a),'Coherence')
    ylabel(ax(a),' Coh','fontsize',strfontsize,'fontname','time new roman')
    xlabel(ax(a),' Hz','fontsize',strfontsize,'fontname','time new roman')
    ax(a).FontSize=strfontsize;
    ax(a).FontName='Time New Roman';


    a=a+1;
    ax(a)=subplot('Position',[.53 .75 .45 .05]);
    title(ax(a),'Acceleration')
    ylabel(ax(a),' g','fontsize',strfontsize,'fontname','time new roman')
    xlabel(ax(a),' seconds','fontsize',strfontsize,'fontname','time new roman')
    ax(a).FontSize=strfontsize;
    ax(a).YLim=.03.*[-1 1];
    ax(a).FontName='Time New Roman';

    a=a+1;
    ax(a)=subplot('Position',[.53 .82 .45 .05]);
    title(ax(a),'FPO7')
    ylabel(ax(a),' Volt','fontsize',strfontsize,'fontname','time new roman')
    ax(a).FontSize=strfontsize;
    ax(a).XTickLabel='';
    ax(a).YLim=5e-4.*[-1 1];
    ax(a).FontName='Time New Roman';

    a=a+1;
    ax(a)=subplot('Position',[.53 .9 .45 .05]);
    title(ax(a),'shear')
    ylabel(ax(a),' Volt','fontsize',strfontsize,'fontname','time new roman')
    ax(a).FontSize=strfontsize;
    ax(a).XTickLabel='';
    ax(a).YLim=.005.*[-1 1];
    ax(a).FontName='Time New Roman';


    [kpan1,Ppan1]=panchev(Profile.epsilon_co(idx1,1),Profile.kvis(idx1));
    [kpan2,Ppan2]=panchev(Profile.epsilon_co(idx1,2),Profile.kvis(idx1));
    [kbatch1,Pbatch1]=batchelor(Profile.epsilon_co(idx1,1), ...
        Profile.chi(idx1), ...
        Profile.kvis(idx1), ...
        .1.*Profile.kvis(idx1));
    [kbatch2,Pbatch2]=batchelor(Profile.epsilon_co(idx1,2), ...
        Profile.chi(idx1,2), ...
        Profile.kvis(idx1), ...
        .1.*Profile.kvis(idx1));
    a=1
    hold(ax(a),'on')
    sca(1)=scatter(ax(a),Profile.epsilon_co(idx1,1),Profile.pr(idx1),200,'pk','filled');
    sca(2)=scatter(ax(a),Profile.epsilon_co(idx1,2),Profile.pr(idx1),200,'dk','filled');
    sca(3)=scatter(ax(a),Profile.epsilon_final(idx1),Profile.pr(idx1),200,'dk','filled');
    sca(1).DisplayName=sprintf('s1:%1.2e',Profile.epsilon_co(idx1,1));
    sca(2).DisplayName=sprintf('s2:%1.2e',Profile.epsilon_co(idx1,2));
    sca(3).DisplayName=sprintf('s:%1.2e',Profile.epsilon_final(idx1));
    hold(ax(a),'off')
    ax(a).XScale='log';

    a=a+1;
    hold(ax(a),'on')
    sca(4)=scatter(ax(a),Profile.chi(idx1,1),Profile.pr(idx1),200,'pk','filled');
    sca(5)=scatter(ax(a),Profile.chi(idx1,2),Profile.pr(idx1),200,'dk','filled');
    sca(4).DisplayName=sprintf('t1:%1.2e',Profile.chi(idx1,1));
    sca(5).DisplayName=sprintf('t2:%1.2e',Profile.chi(idx1,2));
    hold(ax(a),'off')
    ax(a).XScale='log';

    a=a+1;
    hold(ax(a),'on')
    sca(6)=scatter(ax(a),T(idx_ctd),Profile.ctd.P(idx_ctd),200,'pk','filled');
    hold(ax1(2),'on')
    sca(7)=scatter(ax1(2),S(idx_ctd),Profile.ctd.P(idx_ctd),200,'dk','filled');
    hold(ax1(2),'off')
    %     sca(5).DisplayName=sprintf('t1:%1.2e',Profile.chi(idx1,1));
    %     sca(6).DisplayName=sprintf('t2:%1.2e',Profile.chi(idx1,2));
    hold(ax(a),'off')
    %     ax(a).XScale='log';

    a=a+1;
    hold(ax(a),'on')
    sca(8)=scatter(ax(a),Profile.w(idx1),Profile.pr(idx1),200,'pk','filled');
    hold(ax(a),'off')


    a=a+1;
    if isfield(Profile,'vnav')
        hold(ax(a),'on')
        sca(9)=scatter(ax(a),Profile.vnav.pitch(idx_vnav),Profile.vnav.pr(idx_vnav),200,'pk','filled');
        hold(ax(a),'off')
    end

    a=a+1;
    if isfield(Profile,'vnav')
        hold(ax(a),'on')
        sca(10)=scatter(ax(a),Profile.vnav.roll(idx_vnav),Profile.vnav.pr(idx_vnav),200,'pk','filled');
        hold(ax(a),'off')
    end

    a=a+1;
    idx_fc1=find(Profile.f>=Profile.sh_fc(idx1,1),1,'first');
    idx_fc2=find(Profile.f>=Profile.sh_fc(idx1,2),1,'first');
    hold(ax(a),'on')
    loglog(ax(a),Profile.f,Profile.Ps_shear_k.s1(idx1,:),'b--','linewidth',.5);
    loglog(ax(a),Profile.f,Profile.Ps_shear_k.s2(idx1,:),'r--','linewidth',.5);
    
    loglog(ax(a),Profile.f,Profile.Ps_shear_co_k.s1(idx1,:),'b','linewidth',2);
    scatter(ax(a),Profile.f(idx_fc1),Profile.Ps_shear_co_k.s1(idx1,idx_fc1),100,"dy",'filled')
    loglog(ax(a),Profile.f,Profile.Ps_shear_co_k.s2(idx1,:),'r','linewidth',2);
    scatter(ax(a),Profile.f(idx_fc2),Profile.Ps_shear_co_k.s2(idx1,idx_fc2),100,"pk",'filled')

    loglog(ax(a),kpan1.*Profile.w(idx1),Ppan1,'k--');
    loglog(ax(a),kpan2.*Profile.w(idx1),Ppan2,'k--');
    grid(ax(a),'on')
    hold(ax(a),'off')
    ax(a).XScale='log';    
    ax(a).YScale='log';    
    try %
        ax(a).YLim  =[.5.*min(Profile.Ps_shear_k.s1(:),[],'omitnan') max(Profile.Ps_shear_k.s1(:),[],'omitnan')];
    catch % if s1 is only nans
        ax(a).YLim  =[.5.*min(Profile.Ps_shear_k.s2(:),[],'omitnan') max(Profile.Ps_shear_k.s2(:),[],'omitnan')];
    end
    ax(a).XLim  =[Profile.f(2) Profile.f(end)];
    %     title(ax(a),'Shear')
    %     ylabel(ax(a),'s^{-2}/cpm','fontsize',20,'fontname','time new roman')
    %     ax(a).FontSize=20;
    %     ax(a).FontName='Time New Roman';

    a=a+1;
    hold(ax(a),'on')
    idx_fc1=find(Profile.f>=Profile.tg_fc(idx1,1),1,'first');
    idx_fc2=find(Profile.f>=Profile.tg_fc(idx1,2),1,'first');
    loglog(ax(a),Profile.f,Profile.Pt_Tg_k.t1(idx1,:),'linewidth',2);
    scatter(ax(a),Profile.f(idx_fc1),Profile.Pt_Tg_k.t1(idx1,idx_fc1),100,"dy",'filled')
    loglog(ax(a),Profile.f,Profile.Pt_Tg_k.t2(idx1,:),'linewidth',2);
    scatter(ax(a),Profile.f(idx_fc2),Profile.Pt_Tg_k.t2(idx1,idx_fc2),100,"pr",'filled')
    loglog(ax(a),kbatch1.*Profile.w(idx1),Pbatch1,'k--');
    loglog(ax(a),kbatch2.*Profile.w(idx1),Pbatch2,'k--');
    grid(ax(a),'on')
    hold(ax(a),'off')
    ax(a).XScale='log';
    ax(a).YScale='log';
    ax(a).XLim  =[Profile.f(2) Profile.f(end)];
    try
        ax(a).YLim  =[.5.*min(Profile.Pt_Tg_k.t1(idx1,:)) ...
            max(Profile.Pt_Tg_k.t1(idx1,:))];
    catch
        ax(a).YLim  =[.5.*min(Profile.Pt_Tg_k.t2(idx1,:)) ...
            max(Profile.Pt_Tg_k.t2(idx1,:))];
    end


    a=a+1;
    hold(ax(a),'on')
    loglog(ax(a),Profile.f,Profile.Pa_g_f.a1(idx1,:),'linewidth',2);
    loglog(ax(a),Profile.f,Profile.Pa_g_f.a2(idx1,:),'linewidth',2);
    loglog(ax(a),Profile.f,Profile.Pa_g_f.a3(idx1,:),'linewidth',2);
    grid(ax(a),'on')
    hold(ax(a),'off')
    ax(a).XScale='log';
    ax(a).YScale='log';
    ax(a).XLim  =[Profile.f(2) Profile.f(end)];

    a=a+1;
    hold(ax(a),'on')
    semilogx(ax(a),Profile.f,Profile.Cs1a3(idx1,:),'linewidth',2);
    semilogx(ax(a),Profile.f,Profile.Cs2a3(idx1,:),'linewidth',2);
    grid(ax(a),'on')
    hold(ax(a),'off')
    ax(a).XScale='log';
    ax(a).XLim  =[Profile.f(2) Profile.f(end)];
    ax(a).YLim  =[0 1];


    epsi_idx = Profile.ind_range_epsi(idx1,1): ...
        Profile.ind_range_epsi(idx1,2);
    ctd_idx  = Profile.ind_range_ctd(idx1,1): ...
        Profile.ind_range_ctd(idx1,1);
    time_axis=86400.*(Profile.epsi.dnum(epsi_idx)-Profile.epsi.dnum(epsi_idx(1)));

    a=a+1;
    strleg=[];
    hold(ax(a),'on')
    if isfield(Profile.epsi,'a1_g')
        plot(ax(a),time_axis,...
            detrend(Profile.epsi.a1_g(epsi_idx),'constant'),'linewidth',2)
        strleg=[strleg sprintf('%1.2fV',mean(Profile.epsi.a1_g(epsi_idx)))];
    end
    if isfield(Profile.epsi,'a2_g')
        plot(ax(a),time_axis,...
            detrend(Profile.epsi.a2_g(epsi_idx),'constant'),'linewidth',2)

        strleg=[strleg sprintf('%1.2fV',mean(Profile.epsi.a2_g(epsi_idx)))];
    end
    if isfield(Profile.epsi,'a3_g')
        plot(ax(a),time_axis,...
            detrend(Profile.epsi.a3_g(epsi_idx),'constant'),'linewidth',2)

        strleg=[strleg sprintf('%1.2fV',mean(Profile.epsi.a3_g(epsi_idx)))];
    end
    hold(ax(a),'off')
    grid(ax(a),'on')
    leg=legend(ax(a),strleg);
    a=a+1;
    hold(ax(a),'on')
    try
        plot(ax(a),time_axis,...
            detrend(Profile.epsi.t1_volt(epsi_idx),'constant'),'linewidth',2)
    catch
        Profile.epsi.t1_volt=Profile.epsi.dnum.*nan;
        plot(ax(a),time_axis,...
            detrend(Profile.epsi.t1_volt(epsi_idx).*nan,'constant'),'linewidth',2)
    end
    try
        plot(ax(a),time_axis,...
            detrend(Profile.epsi.t2_volt(epsi_idx),'constant'),'linewidth',2)
    catch
        Profile.epsi.t2_volt=Profile.epsi.dnum.*nan;
        plot(ax(a),time_axis,...
            detrend(Profile.epsi.t2_volt(epsi_idx),'constant'),'linewidth',2)
    end
    hold(ax(a),'off')
    grid(ax(a),'on')
    legend(ax(a),{sprintf('%1.2fV',mean(Profile.epsi.t1_volt(epsi_idx))),...
        sprintf('%1.2fV',mean(Profile.epsi.t2_volt(epsi_idx)))})

    a=a+1;
    hold(ax(a),'on')
    plot(ax(a),time_axis,...
        detrend(Profile.epsi.s1_volt(epsi_idx),'constant'),'linewidth',2)
    plot(ax(a),time_axis,...
        detrend(Profile.epsi.s2_volt(epsi_idx),'constant'),'linewidth',2)
    hold(ax(a),'off')
    grid(ax(a),'on')
    legend(ax(a),{sprintf('%1.2fV',mean(Profile.epsi.s1_volt(epsi_idx))),...
        sprintf('%1.2fV',mean(Profile.epsi.s2_volt(epsi_idx)))})

    disp("Click on the epsilon profile.")
    figure(figgrid);
    idx=ginput(1);
    select_dnum=abs(GRID.dnum-idx(1));
    id_profile=find(select_dnum==min(select_dnum));
    select_depth=idx(2);

    if (select_depth<GRID.pr(end) &&  ...
        select_depth>GRID.pr(1) &&    ...
        id_profile<length(GRID.dnum) &&    ...
        id_profile>1)

        a=1;
        cla(ax(a+1));
        cla(ax(a+2));
        cla(ax(a+3));
        cla(ax(a+4));
        cla(ax(a+5));
        cla(ax(a+6));
        cla(ax(a+7));
        cla(ax(a+8));
        cla(ax(a+9));
        cla(ax(a+10));
        cla(ax(a+11));
        cla(ax(a+12));
        delete(sca);
        delete(la1)
        delete(sca1)

    else
        plot_flag=0;
    end

end

% end

