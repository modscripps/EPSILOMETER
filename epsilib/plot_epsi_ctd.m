function [] = plot_epsi_ctd(obj,saveFig)

            if nargin<2
                saveFig=0;
            end

%% Set up axes
figure('units','inches','position',[0 0 10 13])

gap = [0.025 0.025];
margV = [0.08 0.05];
margH = [0.1 0.05];

ax(1)=subtightplot(7,1,1,gap,margV,margH);
ax(2)=subtightplot(7,1,2,gap,margV,margH);
ax(3)=subtightplot(7,1,3,gap,margV,margH);
ax(4)=subtightplot(7,1,4,gap,margV,margH);
ax(5)=subtightplot(7,1,5,gap,margV,margH);
ax(6)=subtightplot(7,1,6,gap,margV,margH);
ax(7)=subtightplot(7,1,7,gap,margV,margH);

cols = obj.plot_properties.Colors; 

%% EPSI

plot(ax(1),obj.epsi.epsitime,obj.epsi.t1_volt,'Color',cols.t1,'LineWidth',obj.plot_properties.LineWidth)
hold(ax(1),'on')
plot(ax(1),obj.epsi.epsitime,obj.epsi.t2_volt,'Color',cols.t2,'LineWidth',obj.plot_properties.LineWidth)
hold(ax(1),'on')

plot(ax(2),obj.epsi.epsitime,obj.epsi.s1_volt,'Color',cols.s1,'LineWidth',obj.plot_properties.LineWidth)
hold(ax(2),'on')
plot(ax(2),obj.epsi.epsitime,obj.epsi.s2_volt,'Color',cols.s2,'LineWidth',obj.plot_properties.LineWidth)
hold(ax(2),'on')

plot(ax(3),obj.epsi.epsitime,obj.epsi.a1_g,'Color',cols.a1,'LineWidth',obj.plot_properties.LineWidth)
hold(ax(3),'on')
plot(ax(3),obj.epsi.epsitime,obj.epsi.a2_g,'Color',cols.a2,'LineWidth',obj.plot_properties.LineWidth)
plot(ax(3),obj.epsi.epsitime,obj.epsi.a3_g,'Color',cols.a3,'LineWidth',obj.plot_properties.LineWidth)
hold(ax(3),'on')
legend(ax(1),{'t1','t2'})
legend(ax(2),{'s1','s2'})
legend(ax(3),{'a1','a2','a3'})

ax(1).XTickLabel='';
ax(2).XTickLabel='';
ax(3).XTickLabel='';
title(ax(1),sprintf('%s-%s-%s',strrep(obj.Meta_Data.mission,'_','\_'),...
    strrep(obj.Meta_Data.vehicle_name,'_','\_'),...
    strrep(obj.Meta_Data.deployment,'_','\_')));
ylabel(ax(1),'FPO7 [Volt]','Fontsize',obj.plot_properties.FontSize,'FontName',obj.plot_properties.FontName)
ylabel(ax(2),'Shear [Volt]','Fontsize',obj.plot_properties.FontSize,'FontName',obj.plot_properties.FontName)
ylabel(ax(3),'Accel [g]','Fontsize',obj.plot_properties.FontSize,'FontName',obj.plot_properties.FontName)
%ax(3).XLabel.String = 'seconds';
set(ax(1),'Fontsize',obj.plot_properties.FontSize,'FontName',obj.plot_properties.FontName)
set(ax(2),'Fontsize',obj.plot_properties.FontSize,'FontName',obj.plot_properties.FontName)
set(ax(3),'Fontsize',obj.plot_properties.FontSize,'FontName',obj.plot_properties.FontName)



%% CTD

plot(ax(4),obj.ctd.ctdtime,obj.ctd.P,'Color',cols.P,'LineWidth',obj.plot_properties.LineWidth)
plot(ax(5),obj.ctd.ctdtime,obj.ctd.T,'Color',cols.T,'LineWidth',obj.plot_properties.LineWidth)
plot(ax(6),obj.ctd.ctdtime,obj.ctd.S,'Color',cols.S,'LineWidth',obj.plot_properties.LineWidth)
ax(4).XTickLabel='';
ax(5).XTickLabel='';
ax(6).XTickLabel='';
ylabel(ax(4),'P [dbar]','Fontsize',obj.plot_properties.FontSize,'FontName',obj.plot_properties.FontName)
ylabel(ax(5),'T [˚C]','Fontsize',obj.plot_properties.FontSize,'FontName',obj.plot_properties.FontName)
ylabel(ax(6),'S [psu]','Fontsize',obj.plot_properties.FontSize,'FontName',obj.plot_properties.FontName)
set(ax(4),'Fontsize',obj.plot_properties.FontSize,'FontName',obj.plot_properties.FontName,'YDir','reverse')
set(ax(5),'Fontsize',obj.plot_properties.FontSize,'FontName',obj.plot_properties.FontName)
set(ax(6),'Fontsize',obj.plot_properties.FontSize,'FontName',obj.plot_properties.FontName)

%% FALL SPEED

plot(ax(7),obj.ctd.ctdtime,movmean(obj.ctd.dPdt,100),'Color',cols.dPdt)
ax(7).XLabel.String = 'ctdtime (sec)';
ylabel(ax(7),'dPdt')
ax(7).YLim = [-0.1 max([1,max(movmean(obj.ctd.dPdt,100))])];

%% AXES

linkaxes(ax,'x')
ax(1).XLim = [obj.epsi.epsitime(1)-15,obj.epsi.epsitime(end)+15];
[ax(:).XGrid] = deal('on');

if saveFig
img = getframe(gcf);
imwrite(img.cdata,fullfile(obj.Meta_Data.datapath,'figs/epsi_ctd_timeseries.png'));
savefig(fullfile(obj.Meta_Data.datapath,'figs/epsi_ctd_timeseries.fig'));
end
