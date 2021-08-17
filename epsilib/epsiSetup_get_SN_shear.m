function [Meta_Data] = set_SN_shear(Meta_Data)
% Set shear probe number and get calibration data

fprintf('**** s1 SN (currently SN = %s, cal = %3.2f)',Meta_Data.AFE.s1.SN,Meta_Data.AFE.s1.cal)
Meta_Data.AFE.s1.SN = input(': ','s');
fprintf('**** s2 SN (currently SN = %s, cal = %3.2f)',Meta_Data.AFE.s2.SN,Meta_Data.AFE.s2.cal)
Meta_Data.AFE.s2.SN = input(': ','s');

Meta_Data=mod_som_get_shear_probe_calibration_v2(Meta_Data);

fprintf('**** s1 SN = %s, cal = %3.2f ****\n',Meta_Data.AFE.s1.SN,Meta_Data.AFE.s1.cal)
fprintf('**** s2 SN = %s, cal = %3.2f ****\n',Meta_Data.AFE.s2.SN,Meta_Data.AFE.s2.cal)

save(fullfile(Meta_Data.datapath,'Meta_Data.mat'),'Meta_Data');
