function mod_epsi_plot_raw(Meta_Data)



load(fullfile(Meta_Data.Epsipath,['epsi_' Meta_Data.deployment '.mat'])) 
load(fullfile(Meta_Data.CTDpath,['epsi_' Meta_Data.deployment '.mat']))

ax(1)=subplot(411);
ax(2)=subplot(412);
ax(3)=subplot(413);
ax(4)=subplot(414);


plot(ax(1),aux1time,P)
plot(ax(2),epsitime,t1)
hold(ax(2),'on')
plot(ax(2),epsitime,t2)
hold(ax(2),'off')
legend(ax(2),'t1','t2')

plot(ax(3),epsitime,s1)
hold(ax(3),'on')
plot(ax(3),epsitime,s2)
hold(ax(3),'off')
title(ax(3),'shear')
legend(ax(3),'s1','s2')

plot(ax(4),epsitime,a1)
hold(ax(4),'on')
plot(ax(4),epsitime,a2)
plot(ax(4),epsitime,a3)
hold(ax(4),'off')
title(ax(4),'accell')
legend(ax(4),'a1','a2','a3')
linkaxes(ax,'x')

switch Meta_Data.PROCESS.recording_mod
    case 'SD'
        print('-dpng',fullfile(Meta_Data.SDRAWpath,[Meta_Data.deployment '_raw.png']))
    case 'STREAMING'
        print('-dpng',fullfile(Meta_Data.RAWpath,[Meta_Data.deployment '_raw.png']))
end

