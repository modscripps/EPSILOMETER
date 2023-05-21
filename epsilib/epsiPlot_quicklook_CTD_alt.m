fctd_mat_dir = '/Users/Shared/EPSI_PROCESSING/Processed/0521_fctd_d4/fctd_mat';

ec = epsi_class(fileparts(fctd_mat_dir));

%fctd = concatenate_fctd_data(fctd_mat_dir);
P = ec.f_get_var_timeseries('P');
dPdt = ec.f_get_var_timeseries('dPdt');
dst = ec.f_get_var_timeseries('dst');
T = ec.f_get_var_timeseries('T');
S = ec.f_get_var_timeseries('S');
acc = ec.f_get_var_timeseries('acceleration');
%%
figure('position',[0 0 704  1013])
gap = [0.05 0.02];
marg_v = [0.04 0.04];
marg_h = [0.07 0.05];

ax(1) = subtightplot(6,1,1,gap,marg_v,marg_h);
plot(P.dnum,P.data,'color',ec.plot_properties.Colors.P); 
title('P')
set(gca,'ydir','reverse')
grid on

ax(2) = subtightplot(6,1,2,gap,marg_v,marg_h);
plot(dPdt.dnum,movmean(dPdt.data,100),'.','color',ec.plot_properties.Colors.dPdt); 
title('dPdt - Fall speed')
grid on

ax(3) = subtightplot(6,1,3,gap,marg_v,marg_h);
plot(dst.dnum,dst.data,'.','color',ec.plot_properties.Colors.alt); 
title('dst - Altimeter height (m)')
grid on

ax(4) = subtightplot(6,1,4,gap,marg_v,marg_h);
plot(T.dnum,T.data,'color',ec.plot_properties.Colors.T); 
title('T')
grid on

ax(5) = subtightplot(6,1,5,gap,marg_v,marg_h);
% plot(S.dnum,S.data,'color',ec.plot_properties.Colors.S); 
% ylim([34.8 36.8])
% title('S')
plot(acc.dnum,acc.data(:,1),'color',ec.plot_properties.Colors.a1)
hold on
plot(acc.dnum,acc.data(:,2),'color',ec.plot_properties.Colors.a2)
plot(acc.dnum,acc.data(:,3),'color',ec.plot_properties.Colors.a3)
title('acceleration')
grid on

ax(6) = subtightplot(6,1,6,gap,marg_v,marg_h);
%plot(fctd.time,fctd.uConductivity(:,1),'.','color',ec.plot_properties.Colors.S)
title('microconductivity')
grid on

lp = linkprop(ax(:),'xlim');
for iAx=1:length(ax)
   datetick(ax(iAx),'x','HH:MM','keeplimits'); 
end