ec = epsi_class;

z = ec.f_get_var_timeseries('z');
dzdt = ec.f_get_var_timeseries('dzdt');

load(fullfile(ec.Meta_Data.paths.mat_data,'PressureTimeseries.mat'));

PT = PressureTimeseries;

figure
ax(1) = subtightplot(2,1,1);
plot(z.dnum,z.data);
hold on
ax(1).YDir = 'reverse';

ax(2) = subtightplot(2,1,2);
plot(dzdt.dnum,movmean(dzdt.data,16*5));
hold on
ax(2).YLim = [-4 4];

lp = linkprop([ax(:)],'xlim');

for iAx=1:2
    axes(ax(iAx))
    for ii=1:length(PT.startprof        s1 = PT.dnum(PT.startprof(ii));
        s2 = PT.dnum(PT.endprof(ii));
        f = fill([s1,s2,s2,s1,s1],[-100,-100,3000,3000,-100],'k');
        f.FaceAlpha = 0.1;
        f.EdgeColor = 'none';
    end
    datetick(ax(iAx),'x','keeplimits');
    grid on
end

ax(1).XLim(1) = 7.390261628205129e+05;
ax(1).YLim(1) = 0;