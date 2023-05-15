ec = epsi_class;
P = ec.f_get_var_timeseries('P');
dPdt = ec.f_get_var_timeseries('dPdt');
dst = ec.f_get_var_timeseries('dst');
T = ec.f_get_var_timeseries('T');
S = ec.f_get_var_timeseries('S');


figure('position',[0 0 704  1013])
ax(1) = subtightplot(5,1,1);
plot(P.dnum,P.data,'color',ec.plot_properties.Colors.P); 
title('P')
set(gca,'ydir','reverse')

ax(2) = subtightplot(5,1,2);
plot(dPdt.dnum,movmean(dPdt.data,100),'.','color',ec.plot_properties.Colors.dPdt); 
title('dPdt - Fall speed')

ax(3) = subtightplot(5,1,3);
plot(dst.dnum,dst.data,'.','color',ec.plot_properties.Colors.alt); 
title('dst - Altimeter height (m)')

ax(4) = subtightplot(5,1,4);
plot(T.dnum,T.data,'color',ec.plot_properties.Colors.T); 
title('T')

ax(5) = subtightplot(5,1,5);
plot(S.dnum,S.data,'color',ec.plot_properties.Colors.S); 
ylim([30,35])
title('S')

lp = linkprop(ax(:),'xlim');
for iAx=1:length(ax)
   datetick(ax(iAx),'x','HH:MM:SS','keeplimits'); 
end