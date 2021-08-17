function plot_profile(Meta_Data,MS,EpsiProfile,CTDProfile)

channels=Meta_Data.PROCESS.channels;
fontsize=20;
figure('units','inch','position',[50,0,15,25]);
% Horizontal acceleration a1 a3
ax(1)=subplot('Position',[.3 .83 .6 .1]);
% Vertical acceleration a2
ax(2)=subplot('Position',[.3 .72 .6 .1]);
% shear s1 s2
ax(3)=subplot('Position',[.3 .61 .6 .1]);
% FPO7 t1 t2
ax(4)=subplot('Position',[.3 .5 .6 .1]);
% spectra
ax(5)=subplot('Position',[.3 .05 .6 .4]);
% CTD Temperature and Salinty profil
ax(6)=subplot('Position',[.1 .1 .1 .8]);


etime=(EpsiProfile.epsitime-EpsiProfile.epsitime(1))*86400;
Dt=etime(end)-etime(1);
tplot=[0 0];splot=[0 0];aplot=[0 0];a2plot=[0 0];
cmap=colormap(lines(length(channels)));
tind=1;aind=1;shind=1;
for c=1:length(channels)
    wh_channels=channels{c};
    switch wh_channels
        case {'t1','t2'}
            hold(ax(4),'on')
            tplot(tind)=plot(ax(4),etime,EpsiProfile.(wh_channels),'Color',cmap(c,:));
            hold(ax(4),'off')
            ax(4).XTick=0:Dt/4:Dt;
            ax(4).XLim=[0 Dt];
            tind=tind+1;
        case {'s1','s2'}
            hold(ax(3),'on')
            splot(shind)=plot(ax(3),etime,EpsiProfile.(wh_channels),'Color',cmap(c,:));
            hold(ax(3),'off')
            ax(3).XTick=0:Dt/4:Dt;
            ax(3).XTickLabel='';
            ax(3).XLim=[0 Dt];
            shind=shind+1;
        case {'a1','a3'}
            hold(ax(1),'on')
            aplot(aind)=plot(ax(1),etime,EpsiProfile.(wh_channels),'Color',cmap(c,:));
            hold(ax(1),'off')
            ax(1).XTick=0:Dt/4:Dt;
            ax(1).XTick=0:Dt/4:Dt;
            aind=aind+1;
        case {'a2'}
            hold(ax(2),'on')
            a2plot=plot(ax(2),etime,EpsiProfile.(wh_channels),'Color',cmap(c,:));
            hold(ax(2),'off')
            ax(2).XTick=0:Dt/4:Dt;
            ax(2).XTickLabel='';
            ax(2).XLim=[0 Dt];
    end
    hold(ax(5),'on')
    loglog(ax(5),MS.f,median(squeeze(MS.Pf(c,:,:))),'Color',cmap(c,:))
    hold(ax(5),'off')
end
title(ax(1),Meta_Data.mission,'Fontsize',fontsize);
xlabel(ax(4),'sec','Fontsize',fontsize)
xlabel(ax(5),'Hz','fontsize',fontsize)
ylabel(ax(5),'s^{-2} /Hz','fontsize',fontsize)
ax(5).XScale='Log';ax(5).YScale='Log';
legend(ax(1),aplot,{'a1','a3'});
legend(ax(2),a2plot,{'a2'});
legend(ax(3),splot,{'s1','s2'});
legend(ax(4),tplot,{'t1','t2'});
grid(ax(5),'on')
ax(1).FontSize=fontsize-2;
ax(2).FontSize=fontsize-2;
ax(3).FontSize=fontsize-2;
ax(4).FontSize=fontsize-2;
ax(5).FontSize=fontsize-2;
ylabel(ax(1),'g','fontsize',fontsize)
ylabel(ax(2),'g','fontsize',fontsize)
ylabel(ax(3),'Volt','fontsize',fontsize)
ylabel(ax(4),'Volt','fontsize',fontsize)



a=6;
[ax1,hl1,hl2]=plotxx(CTDProfile.T,CTDProfile.P,CTDProfile.S,CTDProfile.P,{'',''},{'',''},ax(a));
hl1.Marker='d';
hl2.Marker='d';
hold(ax1(1),'on')
ax1(1).YDir='reverse';
set(ax1(1),'Xscale','linear','Yscale','linear')
ax1(1).XTickLabelRotation=25;
set(ax1(1),'fontsize',15)
xlabel(ax1(1),'CTD T (C) ','fontsize',fontsize)
ylabel(ax1(1),'Pr (dBar)')
grid(ax1(1),'on')
set(ax1(2),'Xscale','linear','Yscale','linear')
ax1(2).XTickLabelRotation=25;
set(ax1(2),'fontsize',15)
ax1(2).YDir='reverse';

xlabel(ax1(2),'CTD S (psu) ','fontsize',fontsize)
