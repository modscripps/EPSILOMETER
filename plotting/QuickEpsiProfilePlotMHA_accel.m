function ax=QuickEpsiProfilePlotMHA_accel(Profile,params)
%function ax=QuickEpsiProfilePlotMHA_accel(Profile,params)
%MHA 8/2024.
%As QuickEpsiProfilePlotMHA but with accel. Do a 6-panel plot of eps, chi, w, T and S for basic profile health.
%basic plotting params can be specified by passing in a structure of the following form:
%    params.epslims=[1e-10 1e-7];
%    params.chilims=[1e-10 1e-6];
%    params.wlims=[.5 .8];
%    params.alims=[-.05 .05];

if nargin < 2
    params.epslims=[1e-10 1e-7];
    params.chilims=[1e-10 1e-6];
    params.wlims=[.5 .8];
    params.alims=[-.05 .05];
    
end

%params.ylims
clf
ax=MySubplot(.1,0.02,0.02,.1,.15,0,6,1);
axes(ax(1))
semilogx(Profile.epsilon_co(:,1),Profile.z,...
    Profile.epsilon_co(:,2),Profile.z);axis ij
xlim(params.epslims)
grid
legend('s1','s2','location','southeast')
xlabel('\epsilon (W/kg)')
ylabel('depth (m)')
axes(ax(2))
semilogx(Profile.chi,Profile.z);axis ij
xlim(params.chilims)
grid
legend('t1','t2','location','southeast')
ax(2).YTickLabel='';
xlabel('\chi (K^2/s)')
axes(ax(3))
plot(Profile.w,Profile.z);axis ij
xlim([0 1])
xlim(params.wlims)
grid
% title([Profile.deployment ', #' num2str(Profile.profNum) ', ' datestr(nanmin(Profile.dnum))],'Interpreter', 'none')
title(['#' num2str(Profile.profNum) ', ' datestr(nanmin(Profile.dnum))],'Interpreter', 'none')
%title([Profile.deployment ', #' num2str(Profile.profNum)],'Interpreter', 'none')
ax(3).YTickLabel='';
xlabel('w (m/s)')
axes(ax(4))
%Put accel onto the Profile grid
a1=interp1(Profile.epsi.dnum,Profile.epsi.a1_g,Profile.dnum);
a2=interp1(Profile.epsi.dnum,Profile.epsi.a2_g,Profile.dnum);
a3=interp1(Profile.epsi.dnum,Profile.epsi.a3_g,Profile.dnum);
a1m=nanmean(a1);
a2m=nanmean(a2);
a3m=nanmean(a3);

plot(a1-a1m,Profile.z,a2-a2m,Profile.z,a3-a3m,Profile.z);axis ij
lg=legend(['a1-' num2str(a1m,1)],['a2-' num2str(a2m)],['a3-' num2str(a3m)],'location','southeast');
pos=get(lg,'position');pos(2)=pos(2)+.25;set(lg,'position',pos);

ax(4).YTickLabel='';
xlim(params.alims)
grid
xlabel('acceleration (g)')

axes(ax(5))
plot(Profile.ctd.T,Profile.ctd.z);axis ij
%xlim([0 1])
grid
ax(5).YTickLabel='';

xlabel('T (^oC)')
axes(ax(6))
plot(Profile.ctd.S,Profile.ctd.z);axis ij
%xlim([0 1])
grid
ax(6).YTickLabel='';

xlabel('S (psu)')
linkaxes(ax,'y');

shg
setpp(15,6)
