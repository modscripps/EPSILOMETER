figure('units','inches','position',[0    0.8472   11.6667   12.2917])

ax(1) = subtightplot(6,1,1);
plot(ctd.dnum,ctd.P);
ax(1).YDir = 'reverse';
legend('P')

ax(2) = subtightplot(6,1,2);
plot(ctd.dnum,movmean(ctd.dPdt,100));
legend('dPdt')

ax(3) = subtightplot(6,1,3);
plot(epsi.dnum,epsi.s1_volt);
hold on
plot(epsi.dnum,epsi.s2_volt);
legend('s1','s2')

ax(4) = subtightplot(6,1,4);
plot(epsi.dnum,epsi.t1_volt);
hold on
plot(epsi.dnum,epsi.t2_volt);
legend('t1','t2')

ax(5) = subtightplot(6,1,5);
plot(epsi.dnum,epsi.a1_g);
legend('a1')

ax(6) = subtightplot(6,1,6);
plot(epsi.dnum,epsi.a2_g);
hold on
plot(epsi.dnum,epsi.a3_g);
legend('a2','a3')

lp = linkprop([ax(:)],'xlim');
for a=1:6
    datetick(ax(a),'x','HH:MM','keeplimits')
end


%ax(2).YLim = [0.42 0.52];
ax(3).YLim = [-0.04 0.04];
ax(5).YLim = [0.99 1.01];
ax(6).YLim = [-0.1 0.1];