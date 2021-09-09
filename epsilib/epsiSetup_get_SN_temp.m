function [Meta_Data] = set_SN_temp(Meta_Data)
% Set temp probe numbers

fprintf('**** t1 SN (currently SN = %s)',Meta_Data.AFE.t1.SN)
Meta_Data.AFE.t1.SN = input(': ','s');
fprintf('**** t2 SN (currently SN = %s)',Meta_Data.AFE.t2.SN)
Meta_Data.AFE.t2.SN = input(': ','s');

fprintf('**** t1 SN = %s ****\n',Meta_Data.AFE.t1.SN)
fprintf('**** t2 SN = %s ****\n',Meta_Data.AFE.t2.SN)

save(fullfile(Meta_Data.paths.data,'Meta_Data.mat'),'Meta_Data');
