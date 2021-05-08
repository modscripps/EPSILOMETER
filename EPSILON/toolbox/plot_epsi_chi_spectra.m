function [ax1,ax2]=plot_epsi_chi_spectra(Profile,scan,ax1,ax2)

if nargin<3
    ax1=alb_create_subplot(2,2,.4,.1,.04,.05);
    ax2=alb_create_subplot_v2(4,1,.05,.1,.38,1,.02,.05);
else
    for i=1:length(ax1)
        hold(ax1(i),'on');
    end
    for i=1:length(ax2)
        hold(ax2(i),'on');
    end
end

%%
G=9.81;
Fn    = .5*325;  % Nyquist frequency
FR    = 2.5;     % Full range in Volts
def_noise=@(x)((FR/2^x)^2 /Fn);
% def_noise20=G ./(Sv2.*scan.w)*def_noise(20)+0*scan.fe;
% def_noise16=def_noise(16)+0*scan.fe;

Accelnoise=G^2*45e-6^2+0*scan.fe;
% Accelnoise_mmp=35e-6^2+0*MS{wh}.f;

%%

title(ax1(3),'TG')
title(ax1(1),'shear')
title(ax1(4),'Accel')
title(ax1(2),'Coh')


a=3;
loglog(ax1(a),scan.fe,scan.Pt1,'b','linewidth',2)
hold(ax1(a),'on')
loglog(ax1(a),scan.fe,scan.Pt2,'r','linewidth',2)
hold(ax1(a),'off')
legend(ax1(a),'t1','t2')
grid(ax1(a),'on')
ax1(a).XLim=[.3 160];
ax1(a).YLim=[1e-10 1e-4];
ax1(a).FontSize=20;
ax1(a).FontName='Times New Roman';
str_title=sprintf('Pr:%3.1f db, N:%3u, NFFT:%3u, $epsilon:%1.2e W kg^{1}$, $chi:%1.2e K^2 s^{-1}$',...
                    scan.pr,scan.N,scan.nfft,scan.epsilon1,scan.chi1);
splt_str_title=strsplit(str_title,',');
splt_str_title{4}=[splt_str_title{4}(1:2) '\' splt_str_title{4}(3:end)];
splt_str_title{5}=[splt_str_title{5}(1:2) '\' splt_str_title{5}(3:end)];
str_title=[splt_str_title{:}];
text(1,2e-4,str_title,'Parent',ax1(a),'fontsize',20,'interpreter','latex')
ylabel(ax1(a),'C^2 /Hz','FontName','Times New Roman','fontsize',20)


a=1;
loglog(ax1(a),scan.fe,scan.Pv1,'b','linewidth',2)
hold(ax1(a),'on')
loglog(ax1(a),scan.fe,scan.Pv2,'r','linewidth',2)
% loglog(ax1(a),scan.fe,def_noise16,'k--','linewidth',1)
% loglog(ax1(a),scan.fe,def_noise20,'k--','linewidth',1)

hold(ax1(a),'off')
legend(ax1(a),'s1','s2')
grid(ax1(a),'on')
ax1(a).XLim=[.3 160];
ax1(a).YLim=[1e-11 1e-6];
ax1(a).FontSize=20;
ax1(a).FontName='Times New Roman';
ylabel(ax1(a),'m^2 s^{-2} /Hz','FontName','Times New Roman','fontsize',20)


a=4;
loglog(ax1(a),scan.fe,scan.Pa2,'b','linewidth',2)
hold(ax1(a),'on')
loglog(ax1(a),scan.fe,scan.Pa1,'r','linewidth',2)
loglog(ax1(a),scan.fe,scan.Pa3,'k','linewidth',2)
loglog(ax1(a),scan.fe,Accelnoise,'k--','linewidth',1)
hold(ax1(a),'off')
legend(ax1(a),'a2','a1','a3')
grid(ax1(a),'on')
ax1(a).XLim=[.3 160];
ax1(a).YLim=[1e-8 1e-4];
ax1(a).FontSize=20;
ax1(a).FontName='Times New Roman';
xlabel(ax1(a),'\omega (Hz)','FontName','Times New Roman','fontsize',20)
ylabel(ax1(a),'m^2 s^{-4} /Hz','FontName','Times New Roman','fontsize',20)


a=2;
semilogx(ax1(a),scan.fe,scan.Cu1a3,'k','linewidth',2)
hold(ax1(a),'on')
% semilogx(ax1(a),scan.fe,scan.Cu2a3,'r','linewidth',2)
semilogx(ax1(a),scan.fe,scan.Cu1a2,'b','linewidth',2)
hold(ax1(a),'off')
legend(ax1(a),'s1a3','s1a2')
grid(ax1(a),'on')
ax1(a).XLim=[.3 160];
ax1(a).FontSize=20;
ax1(a).FontName='Times New Roman';
xlabel(ax1(a),'\omega (Hz)','FontName','Times New Roman','fontsize',20)
ylabel(ax1(a),'coh','FontName','Times New Roman','fontsize',20)


a=1;
plot(ax2(a),Profile.dPdt,Profile.P,'b','linewidth',2)
hold(ax2(a),'on')
plot(ax2(a),[min(Profile.dPdt) max(Profile.dPdt)],[scan.pr scan.pr],'k--','linewidth',2);
plot(ax2(a),[min(Profile.dPdt) max(Profile.dPdt)],[scan.dP(1) scan.dP(1)],'k-.','linewidth',1);
plot(ax2(a),[min(Profile.dPdt) max(Profile.dPdt)],[scan.dP(2) scan.dP(2)],'k-.','linewidth',1);
hold(ax2(a),'off')
legend(ax2(a),'W')
grid(ax2(a),'on')
ax2(a).YLim=[min(Profile.P) max(Profile.P)];
ax2(a).XLim=[min(Profile.dPdt) max(Profile.dPdt)];
ax2(a).YTickLabels='';
axis(ax2(a),'ij')
ax2(a).FontSize=20;
ax2(a).FontName='Times New Roman';
xlabel(ax2(a),'m s^{-1}','FontName','Times New Roman','fontsize',20)


a=2;
[axTS,hl1,hl2]=plotxx(Profile.T,Profile.P,Profile.S,Profile.P,{'T (C)','S (psu)'},{'', ''},ax2(a));
grid(axTS(1),'on')
axTS(1).YLim=[min(Profile.P) max(Profile.P)];
axTS(1).XLim=[min(Profile.T) max(Profile.T)];
axTS(2).YLim=[min(Profile.P) max(Profile.P)];
axTS(2).XLim=[min(Profile.S) max(Profile.S)];
hold(axTS(1),'on')
plot(axTS(1),[min(Profile.T(:)) max(Profile.T(:))],[scan.pr scan.pr],'k--','linewidth',2);
plot(axTS(1),[min(Profile.T(:)) max(Profile.T(:))],[scan.dP(1) scan.dP(1)],'k-.','linewidth',1);
plot(axTS(1),[min(Profile.T(:)) max(Profile.T(:))],[scan.dP(2) scan.dP(2)],'k-.','linewidth',1);
hold(axTS(1),'off')
hl1.LineWidth=2;
hl2.LineWidth=2;
axTS(1).YTickLabels='';
axTS(2).YTickLabels='';
axis(axTS(1),'ij')
axis(axTS(2),'ij')
axTS(1).FontSize=20;
axTS(1).FontName='Times New Roman';
axTS(2).FontSize=20;
axTS(2).FontName='Times New Roman';


a=3;
hepsi=plot(ax2(a),Profile.epsilon,Profile.pr,'linewidth',2);
hold(ax2(a),'on')
plot(ax2(a),[min(Profile.epsilon(:)) max(Profile.epsilon(:))],[scan.pr scan.pr],'k--','linewidth',2);
plot(ax2(a),[min(Profile.epsilon(:)) max(Profile.epsilon(:))],[scan.dP(1) scan.dP(1)],'k-.','linewidth',1);
plot(ax2(a),[min(Profile.epsilon(:)) max(Profile.epsilon(:))],[scan.dP(2) scan.dP(2)],'k-.','linewidth',1);
hold(ax2(a),'off')
legend(ax2(a),'s1','s2')
grid(ax2(a),'on')
ax2(a).XScale='log';
ax2(a).YLim=[min(Profile.P) max(Profile.P)];
ax2(a).XLim=[min(Profile.epsilon(:)) max(Profile.epsilon(:))];
ax2(a).YTickLabels='';
hepsi(1).Color='b';
hepsi(2).Color='r';
axis(ax2(a),'ij')
ax2(a).FontSize=20;
ax2(a).FontName='Times New Roman';
xlabel(ax2(a),'\epsilon (W kg^{-1})','FontName','Times New Roman','fontsize',20)


a=4;
hchi=plot(ax2(a),Profile.chi,Profile.pr,'linewidth',2);
hold(ax2(a),'on')
plot(ax2(a),[min(Profile.chi(:)) max(Profile.chi(:))],[scan.pr scan.pr],'k--','linewidth',2);
plot(ax2(a),[min(Profile.chi(:)) max(Profile.chi(:))],[scan.dP(1) scan.dP(1)],'k-.','linewidth',1);
plot(ax2(a),[min(Profile.chi(:)) max(Profile.chi(:))],[scan.dP(2) scan.dP(2)],'k-.','linewidth',1);
hold(ax2(a),'off')
legend(hchi,'t1','t2')
grid(ax2(a),'on')

ax2(a).YLim=[min(Profile.P) max(Profile.P)];
ax2(a).XLim=[min(Profile.chi(:)) max(Profile.chi(:))];
ax2(a).XScale='log';
hchi(1).Color='b';
hchi(2).Color='r';
axis(ax2(a),'ij')
ylabel(ax2(a),'Depth (m)','fontsize',20,'fontname','Times New Roman')
ax2(a).FontSize=20;
ax2(a).FontName='Times New Roman';
xlabel(ax2(a),'\chi (K^2 s^{-1})','FontName','Times New Roman','fontsize',20)

if nargin>=3
    for i=1:length(ax1)
        hold(ax1(i),'off');
    end
    for i=1:length(ax2)
        hold(ax2(i),'off');
    end
end



