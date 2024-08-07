function epsiPlot_profile_and_spectra(obj,id_profile)

filename=sprintf('Profile%03i.mat',id_profile);
rootpath=strsplit(pwd,'/');
filepath= ...
    fullfile(obj.Meta_Data.paths.profiles,filename);
if isfile(filepath)
    load(filepath);
else
    warning(['Can not find ' filepath])
    return
end

% load('/Volumes/DataDrive1/BLT3/EPSI_PROCESSING/Reprocess/0810_sta08_d15_ef1_pc2/profiles/Profile070.mat')

%
close all
strfontsize=15;

% function plot_epsi_spectra(Profile)
T=filloutliers(Profile.ctd.T,'linear','movmean',floor(length(Profile.ctd.T)/10));
S=filloutliers(real(Profile.ctd.S),'linear','movmean',floor(length(Profile.ctd.T)/10));
% T=Profile.ctd.T;
% S=Profile.ctd.S;

if isfield(Profile,'vnav')
    Profile.vnav=vnav_pitch_roll_heading(Profile.vnav);
    Profile.vnav.pr=interp1(Profile.ctd.dnum,Profile.ctd.P,Profile.vnav.dnum);
else
    Profile.vnav.pitch   = Profile.pr.* nan;
    Profile.vnav.roll    = Profile.pr.* nan;
    Profile.vnav.heading = Profile.pr.* nan;
    Profile.vnav.dnum    = Profile.dnum;
    Profile.vnav.pr      = Profile.pr;
end

%% Make figure
fig = figure;
% Set figure size based on screen size
defaultFigWidth = 1680;
defaultFigHeight = 886;
screenSize = get(0,'screensize');
mult = round(min([screenSize(3)/defaultFigWidth,screenSize(4)/defaultFigHeight]),2);
set(fig,'Units','pixels','Position',[1 1 defaultFigWidth*mult defaultFigHeight*mult]);

ax(1)=subplot('Position',[0.04 0.08 0.08 0.8]); %epsilon
ax(2)=subplot('Position',[.13 .08 .08 .8]); %chi
ax(3)=subplot('Position',[.22 .08 .05 .8]); %T and S
ax(4)=subplot('Position',[.4 .08 .05 .8]);  %pitch
ax(5)=subplot('Position',[.28 .08 .05 .8]); %roll
ax(6)=subplot('Position',[.34 .08 .05 .8]); %w
ax(7)=subplot('Position',[.53 .5 .45 .13]); %shear spectra
ax(8)=subplot('Position',[0.53 0.3477 0.45 0.13]); %fpo7 spectra
ax(9)=subplot('Position',[.53 .19 .45 .13]); %acceleation spectra
ax(10)=subplot('Position',[.53 .08 .45 .08]); %coherence
ax(11)=subplot('Position',[.53 .7 .45 .07]); %shear
ax(12)=subplot('Position',[.53 .8 .45 .07]); %fpo7
ax(13)=subplot('Position',[.53 .9 .45 .07]); %acceleration

%%

% Axes 1
a=1;
h = semilogx(ax(a),Profile.epsilon_co,Profile.pr);
h(1).Color = obj.plot_properties.Colors.s1;
h(2).Color = obj.plot_properties.Colors.s2;
hold(ax(a),'on')
semilogx(ax(a),Profile.epsilon_final,Profile.pr,'k','linewidth',2);
hold(ax(a),'off')
leg(a)=legend(ax(a),'s1','s2','location','northoutside');
leg(a).Position=[0.0436 0.9016 0.0732 0.0807];
ax(a).FontSize=strfontsize;
ax(a).FontName='Time New Roman';
grid(ax(a),'on')
axis(ax(a),'ij')
ax(a).XLim=[min(Profile.epsilon_co(:)) max(Profile.epsilon_co(:))];
ax(a).YLim=[max(min(Profile.pr),0) max(Profile.pr)];
xlabel(ax(a),'\epsilon','fontsize',strfontsize,'fontname','time new roman')
ylabel(ax(a),'Depth [m]','fontsize',strfontsize,'fontname','time new roman')

% Axes 2
a=2;
h = semilogx(ax(a),Profile.chi,Profile.pr);
h(1).Color = obj.plot_properties.Colors.t1;
h(2).Color = obj.plot_properties.Colors.t2;
leg(a)=legend(ax(a),'t1','t2','location','northoutside');
leg(a).Position=[0.1327 0.9016 0.0723 0.0807];
ax(a).FontSize=strfontsize;
ax(a).FontName='Time New Roman';
grid(ax(a),'on')
axis(ax(a),'ij')
try
ax(a).XLim=[min(Profile.chi(:)) max(Profile.chi(:))];
catch
ax(a).XLim=[-10 -3];
end
ax(a).YLim=[max(min(Profile.pr),0) max(Profile.pr)];
ax(a).YTickLabel='';
xlabel(ax(a),'\chi','fontsize',strfontsize,'fontname','time new roman')

% Axes 3
a=3;
axes(ax(3))
ax1(1) = axes('position',ax(3).Position);
h(1) = plot(T,Profile.ctd.P);
ax1(2) = axes('position',ax(3).Position);
ax1(2).Color = 'none';
h(2) = plot(S,Profile.ctd.P);
ax1(1).XAxisLocation = 'top';
h(1).Color = obj.plot_properties.Colors.T;
h(2).Color = obj.plot_properties.Colors.S;
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
ax1(1).XColor = obj.plot_properties.Colors.T;
ax1(2).XColor = obj.plot_properties.Colors.S;
xlabel(ax(a),'T','fontsize',strfontsize,'fontname','time new roman','color',obj.plot_properties.Colors.T)
xlabel(ax1(2),'S','fontsize',strfontsize,'fontname','time new roman','color',obj.plot_properties.Colors.S)
ax1(2).FontSize=strfontsize;


% Axes 4
a=4;
plot(ax(a),Profile.w,Profile.pr,'color',obj.plot_properties.Colors.dPdt);
axis(ax(a),'ij')
ax(a).FontSize=strfontsize;
ax(a).FontName='Time New Roman';
grid(ax(a),'on')
axis(ax(a),'ij')
ax(a).XLim=[min(Profile.w) max(Profile.w)];
ax(a).YLim=[max(min(Profile.pr),0) max(Profile.pr)];
ax(a).YTickLabel='';
xlabel(ax(a),'w [m/s]','fontsize',strfontsize,'fontname','time new roman')

% Axes 5
a=5;
if isfield(Profile,'vnav')
plot(ax(a),Profile.vnav.pitch,Profile.vnav.pr);
axis(ax(a),'ij')
ax(a).FontSize=strfontsize;
ax(a).FontName='Time New Roman';
grid(ax(a),'on')
axis(ax(a),'ij')
try
    ax(a).XLim=[min(Profile.vnav.pitch) max(Profile.vnav.pitch)];
catch
    ax(a).XLim=[0 1];
end
ax(a).YLim=[max(min(Profile.pr),0) max(Profile.pr)];
ax(a).YTickLabel='';
xlabel(ax(a),'pitch','fontsize',strfontsize,'fontname','time new roman')
end

% Axes 6
a=6;
if isfield(Profile,'vnav')
plot(ax(a),unwrap(Profile.vnav.roll)+pi,Profile.vnav.pr);
axis(ax(a),'ij')
ax(a).FontSize=strfontsize;
ax(a).FontName='Time New Roman';
str_title=sprintf('%s %s Profile %i ',            ...
                  Profile.Meta_Data.mission,      ...
                  Profile.Meta_Data.vehicle_name, ...
                  Profile.profNum);

text(ax(a),-pi/2,Profile.vnav.pr(1)-15,{rootpath{end},...
    str_title},...
    'fontsize',strfontsize,...
    'fontname','time new roman',...
    'Interpreter','none')
grid(ax(a),'on')
axis(ax(a),'ij')
try
ax(a).XLim=[min(unwrap(Profile.vnav.roll)+pi) max(unwrap(Profile.vnav.roll)+pi)];
catch
    ax(a).XLim=[0 1];
end
ax(a).YLim=[max(min(Profile.pr),0) max(Profile.pr)];
ax(a).YTickLabel='';
xlabel(ax(a),'roll','fontsize',strfontsize,'fontname','time new roman')
end

% Axes 7
a=7;
title(ax(a),'Shear','fontsize',strfontsize,'fontname','time new roman')
ylabel(ax(a),'s^{-2}/cpm','fontsize',strfontsize,'fontname','time new roman')
grid(ax(a),'on')
ax(a).FontSize=strfontsize;
ax(a).XTickLabel='';
ax(a).FontName='Time New Roman';

% Axes 8
a=8;
title(ax(a),'FPO7')
ylabel(ax(a),'˚K^2 m^{-2}/cpm','fontsize',strfontsize,'fontname','time new roman')
grid(ax(a),'on')
ax(a).XTickLabel='';
ax(a).FontSize=strfontsize;
ax(a).FontName='Time New Roman';

% Axes 9
a=9;
title(ax(a),'Acceleration')
ylabel(ax(a),' g^2 /cpm','fontsize',strfontsize,'fontname','time new roman')
ax(a).XTickLabel='';
grid(ax(a),'on')
ax(a).FontSize=strfontsize;
ax(a).FontName='Time New Roman';

% Axes 10
a=10;
grid(ax(a),'on')
title(ax(a),'Coherence')
ylabel(ax(a),' Coh','fontsize',strfontsize,'fontname','time new roman')
xlabel(ax(a),' Hz','fontsize',strfontsize,'fontname','time new roman')
ax(a).FontSize=strfontsize;
ax(a).FontName='Time New Roman';

% Axes 11
a=11;
title(ax(a),'Acceleration')
ylabel(ax(a),' g','fontsize',strfontsize,'fontname','time new roman')
xlabel(ax(a),' seconds','fontsize',strfontsize,'fontname','time new roman')
ax(a).FontSize=strfontsize;
%ax(a).YLim=.03.*[-1 1];
ax(a).FontName='Time New Roman';

% Axes 12
a=12;
title(ax(a),'FPO7')
ylabel(ax(a),' Volt','fontsize',strfontsize,'fontname','time new roman')
ax(a).FontSize=strfontsize;
ax(a).XTickLabel='';
%ax(a).YLim=5e-4.*[-1 1];
ax(a).FontName='Time New Roman';

% Axes 13
a=13;
title(ax(a),'shear')
ylabel(ax(a),' Volt','fontsize',strfontsize,'fontname','time new roman')
ax(a).FontSize=strfontsize;
ax(a).XTickLabel='';
%ax(a).YLim=.005.*[-1 1];
ax(a).FontName='Time New Roman';

plot_flag=1;
first_time=1;
while plot_flag
    disp("Click on the epsilon profile.")
    if first_time==1   
        idx=ginput(1);
        first_time=0;
    end
    
    idx1=find(Profile.pr>idx(2),1,'first');
    idx_ctd=find(Profile.ctd.P>idx(2),1,'first');
    idx_vnav=find(Profile.vnav.pr>idx(2),1,'first');
    
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
    a=1;
    hold(ax(a),'on')
    sca(1)=scatter(ax(a),Profile.epsilon_co(idx1,1),Profile.pr(idx1),200,'ok','filled','markerfacecolor',obj.plot_properties.Colors.s1);
    sca(2)=scatter(ax(a),Profile.epsilon_co(idx1,2),Profile.pr(idx1),200,'ok','filled','markerfacecolor',obj.plot_properties.Colors.s2);
    
    sca(1).DisplayName=sprintf('s1:%1.2e',Profile.epsilon(idx1,1));
    sca(2).DisplayName=sprintf('s2:%1.2e',Profile.epsilon(idx1,2));
    hold(ax(a),'off')
    ax(a).XScale='log';    
    
    a=2;
    hold(ax(a),'on')
    sca(3)=scatter(ax(a),Profile.chi(idx1,1),Profile.pr(idx1),200,'ok','filled','markerfacecolor',obj.plot_properties.Colors.t1);
    sca(4)=scatter(ax(a),Profile.chi(idx1,2),Profile.pr(idx1),200,'ok','filled','markerfacecolor',obj.plot_properties.Colors.t2);
    sca(3).DisplayName=sprintf('t1:%1.2e',Profile.chi(idx1,1));
    sca(4).DisplayName=sprintf('t2:%1.2e',Profile.chi(idx1,2));
    hold(ax(a),'off')
    ax(a).XScale='log';    

    a=3;
    hold(ax(a),'on')
    sca(5)=scatter(ax(a),T(idx_ctd),Profile.ctd.P(idx_ctd),200,'ok','filled','markerfacecolor',obj.plot_properties.Colors.T);
    hold(ax1(2),'on')
    sca(6)=scatter(ax1(2),S(idx_ctd),Profile.ctd.P(idx_ctd),200,'ok','filled','markerfacecolor',obj.plot_properties.Colors.S);
    hold(ax1(2),'off')
%     sca(5).DisplayName=sprintf('t1:%1.2e',Profile.chi(idx1,1));
%     sca(6).DisplayName=sprintf('t2:%1.2e',Profile.chi(idx1,2));
    hold(ax(a),'off')
%     ax(a).XScale='log';    

    a=4;
    hold(ax(a),'on')
    sca(7)=scatter(ax(a),Profile.w(idx1),Profile.pr(idx1),200,'ok','filled','markerfacecolor',obj.plot_properties.Colors.dPdt);
    hold(ax(a),'off')

    a=5;
    hold(ax(a),'on')
    sca(8)=scatter(ax(a),Profile.vnav.pitch(idx_vnav),Profile.vnav.pr(idx_vnav),200,'ok','filled');
    hold(ax(a),'off')

    a=6;
    hold(ax(a),'on')
    sca(9)=scatter(ax(a),Profile.vnav.roll(idx_vnav),Profile.vnav.pr(idx_vnav),200,'ok','filled');
    hold(ax(a),'off')


    a=7;
    idx_fc1=find(Profile.f>=Profile.sh_fc(idx1,1),1,'first');
    idx_fc2=find(Profile.f>=Profile.sh_fc(idx1,2),1,'first');
    hold(ax(a),'on')
    loglog(ax(a),Profile.f,Profile.Ps_shear_k.s1(idx1,:),'--','linewidth',.5,'color',obj.plot_properties.Colors.s1);
    loglog(ax(a),Profile.f,Profile.Ps_shear_k.s2(idx1,:),'--','linewidth',.5,'color',obj.plot_properties.Colors.s2);
    
    loglog(ax(a),Profile.f,Profile.Ps_shear_co_k.s1(idx1,:),'linewidth',2,'color',obj.plot_properties.Colors.s1);
    scatter(ax(a),Profile.f(idx_fc1),Profile.Ps_shear_co_k.s1(idx1,idx_fc1),100,"ok",'filled')
    loglog(ax(a),Profile.f,Profile.Ps_shear_co_k.s2(idx1,:),'linewidth',2,'color',obj.plot_properties.Colors.s2);
    scatter(ax(a),Profile.f(idx_fc2),Profile.Ps_shear_co_k.s2(idx1,idx_fc2),100,"ok",'filled')

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

    a=8;
    hold(ax(a),'on')
    idx_fc1=find(Profile.f>=Profile.tg_fc(idx1,1),1,'first');
    idx_fc2=find(Profile.f>=Profile.tg_fc(idx1,2),1,'first');
    loglog(ax(a),Profile.f,Profile.Pt_Tg_k.t1(idx1,:),'linewidth',2,'color',obj.plot_properties.Colors.t1);
    scatter(ax(a),Profile.f(idx_fc1),Profile.Pt_Tg_k.t1(idx1,idx_fc1),100,"dy",'filled')
    loglog(ax(a),Profile.f,Profile.Pt_Tg_k.t2(idx1,:),'linewidth',2,'color',obj.plot_properties.Colors.t2);
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

    
    a=9;
    hold(ax(a),'on')
    loglog(ax(a),Profile.f,Profile.Pa_g_f.a1(idx1,:),'linewidth',2,'color',obj.plot_properties.Colors.a1);
    loglog(ax(a),Profile.f,Profile.Pa_g_f.a2(idx1,:),'linewidth',2,'color',obj.plot_properties.Colors.a2);
    loglog(ax(a),Profile.f,Profile.Pa_g_f.a3(idx1,:),'linewidth',2,'color',obj.plot_properties.Colors.a3);
    grid(ax(a),'on')
    hold(ax(a),'off')
    ax(a).XScale='log';    
    ax(a).YScale='log';    
    ax(a).XLim  =[Profile.f(2) Profile.f(end)];    

    a=10;
    hold(ax(a),'on')
    semilogx(ax(a),Profile.f,Profile.Cs1a3(idx1,:),'linewidth',2,'color',obj.plot_properties.Colors.s1);
    semilogx(ax(a),Profile.f,Profile.Cs2a3(idx1,:),'linewidth',2,'color',obj.plot_properties.Colors.s2);
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

    a=11;
    strleg=[];
    hold(ax(a),'on')
    if isfield(Profile.epsi,'a1_g')
        plot(ax(a),time_axis,...
            detrend(Profile.epsi.a1_g(epsi_idx),'constant'),'linewidth',2,'color',obj.plot_properties.Colors.a1)
        strleg=[strleg sprintf('%1.2fV',mean(Profile.epsi.a1_g(epsi_idx)))];
    end
    if isfield(Profile.epsi,'a2_g')
    plot(ax(a),time_axis,...
        detrend(Profile.epsi.a2_g(epsi_idx),'constant'),'linewidth',2,'color',obj.plot_properties.Colors.a2)

        strleg=[strleg sprintf('%1.2fV',mean(Profile.epsi.a2_g(epsi_idx)))];
    end
    if isfield(Profile.epsi,'a3_g')
    plot(ax(a),time_axis,...
        detrend(Profile.epsi.a3_g(epsi_idx),'constant'),'linewidth',2,'color',obj.plot_properties.Colors.a3)

        strleg=[strleg sprintf('%1.2fV',mean(Profile.epsi.a3_g(epsi_idx)))];
    end
    hold(ax(a),'off')
    grid(ax(a),'on')
    leg=legend(ax(a),strleg);

    a=12;
    hold(ax(a),'on')
    try
        plot(ax(a),time_axis,...
            detrend(Profile.epsi.t1_volt(epsi_idx),'constant'),'linewidth',2,'color',obj.plot_properties.Colors.t1)
    catch
        Profile.epsi.t1_volt=Profile.epsi.dnum.*nan;
        plot(ax(a),time_axis,...
            detrend(Profile.epsi.t1_volt(epsi_idx).*nan,'constant'),'linewidth',2,'color',obj.plot_properties.Colors.t1)
    end
    try
    plot(ax(a),time_axis,...
        detrend(Profile.epsi.t2_volt(epsi_idx),'constant'),'linewidth',2,'color',obj.plot_properties.Colors.t2)
    catch
    Profile.epsi.t2_volt=Profile.epsi.dnum.*nan;
    plot(ax(a),time_axis,...
        detrend(Profile.epsi.t2_volt(epsi_idx),'constant'),'linewidth',2,'color',obj.plot_properties.Colors.t2)
    end
    hold(ax(a),'off')
    grid(ax(a),'on')
    legend(ax(a),{sprintf('%1.2fV',mean(Profile.epsi.t1_volt(epsi_idx))),...
                  sprintf('%1.2fV',mean(Profile.epsi.t2_volt(epsi_idx)))})

    a=13;
    hold(ax(a),'on')
    plot(ax(a),time_axis,...
        detrend(Profile.epsi.s1_volt(epsi_idx),'constant'),'linewidth',2,'color',obj.plot_properties.Colors.s1)
    plot(ax(a),time_axis,...
        detrend(Profile.epsi.s2_volt(epsi_idx),'constant'),'linewidth',2,'color',obj.plot_properties.Colors.s2)
    hold(ax(a),'off')
    grid(ax(a),'on')
    legend(ax(a),{sprintf('%1.2fV',mean(Profile.epsi.s1_volt(epsi_idx))),...
                  sprintf('%1.2fV',mean(Profile.epsi.s2_volt(epsi_idx)))})

    disp("Click on the epsilon profile.")
    idx=ginput(1);
    if (idx(2)<Profile.pr(end) &&  ...
        idx(2)>Profile.pr(1) &&    ...
        idx(1)<max(Profile.epsilon_co(:)) &&    ...
        idx(1)>min(Profile.epsilon_co(:)))
        
        cla(ax(7));
        cla(ax(8));
        cla(ax(9));
        cla(ax(10));
        cla(ax(11));
        cla(ax(12));
        cla(ax(13));
        delete(sca);

    else
        plot_flag=0;    
    end

end

% end


