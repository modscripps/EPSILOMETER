if isfield(epsi,'efe_time_laptop')
    
fig = figure('position',[0 0 10 13]);
gap = [0.08 0.02];
margV = [0.07 0.05];
margH = [0.14 0.05];

ax(1) = subtightplot(3,1,1,gap,margV,margH);
x = find(~epsi.efe_bad_blocks);
colors = ax(1).ColorOrder;

plot(x,epsi.efe_time_sendout(x),'co','Color',colors(1,:))
title('epsi.efe\_time\_sendout');
% xlabel('sample #')
ylabel('time value')

ax(2) = subtightplot(3,1,2,gap,margV,margH);
plot(x,epsi.efe_time_laptop(x),'go')
title('epsi.efe\_time\_laptop');
xlabel('sample #')
ylabel('time value')

ax(3) = subtightplot(3,1,3,gap.*2,margV,margH);
plot(epsi.epsitime,'k*')
title('epsi.epsitime')
xlabel('sample #')
ylabel('time value')

figureStamp(getFilename)
img = getframe(gcf);
imwrite(img.cdata,fullfile(Meta_Data.datapath,['figs/raw_efe_timestamps.png']));


%%
figure
x = 1:length(ctd.timestamp);

[ax,h(1),h(2)] = plotyy(x,ctd.ctdtime,x,ctd.P);
[h(:).Marker] = deal('.');
[h(:).LineStyle] = deal('none');
ax(1).YLabel.String = 'ctd time';
ax(2).YLabel.String = 'ctd prressure';
ax(2).YDir = 'reverse';

ax(1).Position(1) = 0.1;
ax(2).Position(1) = 0.1;
ax(1).Position(3) = 0.75;
ax(2).Position(3) = 0.75;

figureStamp(getFilename)
img = getframe(gcf);
imwrite(img.cdata,fullfile(Meta_Data.datapath,['figs/raw_ctd_timestamps.png']));

end