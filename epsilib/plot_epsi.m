function [ax] = plot_epsi(obj)
% Epsi class function [] = plot_epsi(obj)
%
% INPUTS:
%   obj = epsi_class or structure which must include 'epsi' field
%

if ~any([isfield(obj,'plot_properties'),isclassfield(obj,'plot_properties')])
    obj.plot_properties = epsiSetup_set_plot_properties;
end
if ~any([isfield(obj,'Meta_Data'),isclassfield(obj,'Meta_Data')])
    obj.Meta_Data.mission = '';
    obj.Meta_Data.vehicle_name = '';
    obj.Meta_Data.deployment = '';
end

figure('units','inches','position',[0 0 9.5 12])
ax(1)=subplot(311);
ax(2)=subplot(312);
ax(3)=subplot(313);

% Reference time in seconds to start of record
time = obj.epsi.time_s-obj.epsi.time_s(1);

plot(ax(1),time,obj.epsi.t1_volt,'color',obj.plot_properties.Colors.t1,'LineWidth',obj.plot_properties.LineWidth)
hold(ax(1),'on')
plot(ax(1),time,obj.epsi.t2_volt,'color',obj.plot_properties.Colors.t2,'LineWidth',obj.plot_properties.LineWidth)

plot(ax(2),time,obj.epsi.s1_volt,'color',obj.plot_properties.Colors.s1,'LineWidth',obj.plot_properties.LineWidth)
hold(ax(2),'on')
plot(ax(2),time,obj.epsi.s2_volt,'color',obj.plot_properties.Colors.s2,'LineWidth',obj.plot_properties.LineWidth)

plot(ax(3),time,obj.epsi.a1_g,'color',obj.plot_properties.Colors.a1,'LineWidth',obj.plot_properties.LineWidth)
hold(ax(3),'on')
plot(ax(3),time,obj.epsi.a2_g,'color',obj.plot_properties.Colors.a2,'LineWidth',obj.plot_properties.LineWidth)
plot(ax(3),time,obj.epsi.a3_g,'color',obj.plot_properties.Colors.a3,'LineWidth',obj.plot_properties.LineWidth)

%ALB mean in legend
legend(ax(1),{sprintf('t1 %3.2f V',nanmean(obj.epsi.t1_volt)),...
    sprintf('t2 %3.2f V',nanmean(obj.epsi.t2_volt))})
legend(ax(2),{sprintf('s1 %3.2f V',nanmean(obj.epsi.s1_volt)),...
    sprintf('s2 %3.2f V',nanmean(obj.epsi.s2_volt))})
legend(ax(3),{sprintf('a1 %3.2f g',nanmean(obj.epsi.a1_g)),...
    sprintf('a2 %3.2f g',nanmean(obj.epsi.a2_g)),...
    sprintf('a3 %3.2f g',nanmean(obj.epsi.a3_g))})

ax(1).XTickLabel='';
ax(2).XTickLabel='';
title(ax(1),sprintf('%s-%s-%s   :',strrep(obj.Meta_Data.mission,'_','\_'),...
    strrep(obj.Meta_Data.vehicle_name,'_','\_'),...
    strrep(obj.Meta_Data.deployment,'_','\_')));
ylabel(ax(1),'FPO7 [Volt]','Fontsize',obj.plot_properties.FontSize,'FontName',obj.plot_properties.FontName)
ylabel(ax(2),'Shear [Volt]','Fontsize',obj.plot_properties.FontSize,'FontName',obj.plot_properties.FontName)
ylabel(ax(3),'Accel [g]','Fontsize',obj.plot_properties.FontSize,'FontName',obj.plot_properties.FontName)
ax(3).XLabel.String = 'seconds';
set(ax(1),'Fontsize',obj.plot_properties.FontSize,'FontName',obj.plot_properties.FontName)
set(ax(2),'Fontsize',obj.plot_properties.FontSize,'FontName',obj.plot_properties.FontName)
set(ax(3),'Fontsize',obj.plot_properties.FontSize,'FontName',obj.plot_properties.FontName)
linkaxes(ax,'x')