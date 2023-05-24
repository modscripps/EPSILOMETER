%% Plot temperature, epsilon, and chi with density contours
data=load(fullfile(ec.Meta_Data.paths.profiles,'griddedProfiles'));

data.GRID.bottom_depth=filloutliers(data.GRID.bottom_depth,'linear');
close all
figure('units','inches','position',[0         0   15.3194   13.1111])

dnummask=~isnan(data.GRID.dnum);
ax(1) =subtightplot(3,1,1);

pcolorjw(data.GRID.dnum(dnummask),data.GRID.z,data.GRID.t(:,dnummask));
hold on
%n_fill_bathy(data.GRID.dnum(dnummask),data.GRID.bottom_depth(dnummask))
%[c,ch]=contour(data.GRID.dnum(dnummask),data.GRID.z,data.GRID.sgth(:,dnummask),20,'k');
caxis([9 20])
cax1=colorbar;
grid(ax(1),'on');
%n_fill_bathy(data.GRID.dnum(dnummask),data.GRID.bottom_depth(dnummask))
title('Temperature','fontname','Times New Roman','fontsize',20)
ylabel('Depth','fontname','Times New Roman','fontsize',20)
ylabel(cax1,'Celsius','fontname','Times New Roman','fontsize',20)

ax(2) = subtightplot(3,1,2);
pcolorjw(data.GRID.dnum(dnummask),data.GRID.z,log10(data.GRID.epsilon_final(:,dnummask)));
cax2=colorbar;
title('Epsilon')
ylabel('Depth','fontname','Times New Roman','fontsize',20);
hold on
[c,ch]=contour(data.GRID.dnum(dnummask),data.GRID.z,data.GRID.sgth(:,dnummask),'k','levellist',-10:0.5:-6);
caxis([-10 -6])
grid(ax(2),'on');

%n_fill_bathy(data.GRID.dnum(dnummask),data.GRID.bottom_depth(dnummask))
ylabel(cax2,'W kg^{-1}','fontname','Times New Roman','fontsize',20)

ax(3) = subtightplot(3,1,3);
pcolorjw(data.GRID.dnum(dnummask),data.GRID.z,log10(data.GRID.chi1(:,dnummask)));
hold on
[c,ch]=contour(data.GRID.dnum(dnummask),data.GRID.z,data.GRID.sgth(:,dnummask),'k','levellist',-10:0.5:-5);
%n_fill_bathy(data.GRID.dnum(dnummask),data.GRID.bottom_depth(:,dnummask))
title('Chi 1')
cax3=colorbar;
grid(ax(3),'on');

caxis([-10 -5])
ylabel('Depth','fontname','Times New Roman','fontsize',20)
xlabel(datestr(data.GRID.dnum(1),"dd-mm"),'fontname','Times New Roman','fontsize',20)
[ax(:).YDir] = deal('reverse');
for a=1:length(ax)
    ax(a).XTick=data.GRID.dnum(1):2/24:data.GRID.dnum(end);
    ax(a).XTickLabel='';
end
ax(3).XTickLabel=datestr(data.GRID.dnum(1):2/24:data.GRID.dnum(end),'DD - HH:MM');
ax(3).XTickLabelRotation=45;
ylabel(cax3,'K^2 s^{-1}','fontname','Times New Roman','fontsize',20)

lp = linkprop([ax(:)],'xlim');

% Save section plot
save_name = fullfile(ec.Meta_Data.paths.figures,'deployment_sections');
eval(['export_fig ' save_name ' -png -r150 -nocrop']);
%print('-dpng2',save_name)
