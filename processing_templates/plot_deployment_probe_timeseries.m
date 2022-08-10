% -------------------------------------------
depl_dir = '0808_sta07_d13_ef1_pc2';
% -------------------------------------------

lab_dir = fullfile('/Volumes/FCTD_EPSI/Deployments/',depl_dir);
processing_dir = fullfile('/Users/Shared/EPSI_PROCESSING/',depl_dir);
md = '/Volumes/FCTD Softwares used in BLT 2022/EPSILOMETER_FCTD/Meta_Data_Process/Meta_Data_Process_blt_2022.txt';

% % Rsync deployment files to the processing directory
% com = sprintf('/usr/bin/rsync -auv %s %s',lab_dir,processing_dir);  %
% unix(com);

ec1 = epsi_class(processing_dir,md);
t1 = ec1.f_get_var_timeseries('t1_volt');
t2 = ec1.f_get_var_timeseries('t2_volt');
s1 = ec1.f_get_var_timeseries('s1_volt');
s2 = ec1.f_get_var_timeseries('s2_volt');
P  = ec1.f_get_var_timeseries('P');
T  = ec1.f_get_var_timeseries('T');

figure('position',[610     1   804   984])
margH = [0.05 0.05];
margV = [0.05 0.05];
gap   = [0.025 0.05];

ax(1) = subtightplot(6,1,1,gap,margV,margH);
plot(P.dnum,P.data,'color',ec1.plot_properties.Colors.P);
ylabel('Pressure')
ax(1).YDir = 'reverse';

ax(2) = subtightplot(6,1,2,gap,margV,margH);
plot(s1.dnum,s1.data,'color',ec1.plot_properties.Colors.s1);
ylabel('s1')

ax(3) = subtightplot(6,1,3,gap,margV,margH);
plot(s2.dnum,s2.data,'color',ec1.plot_properties.Colors.s2);
ylabel('s2')

ax(4) = subtightplot(6,1,4,gap,margV,margH);
plot(t1.dnum,t1.data,'color',ec1.plot_properties.Colors.t1);
ylabel('t1')

ax(5) = subtightplot(6,1,5,gap,margV,margH);
plot(t2.dnum,t2.data,'color',ec1.plot_properties.Colors.t2);
ylabel('t2')

ax(6) = subtightplot(6,1,6,gap,margV,margH);
plot(T.dnum,T.data,'color',ec1.plot_properties.Colors.T);
ylabel('CTD T')


for a=1:length(ax)
   datetick(ax(a),'x','HH:MM','keeplimits'); 
end

[ax(1:5).XTickLabel] = deal('');
lp = linkprop([ax(:)],'xlim');