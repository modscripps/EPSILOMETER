function [Meta_Data] = set_SN_temp(Meta_Data)
% Set temp probe numbers

% NC 10/7/21 - Check for 'AFE' or 'epsi' strucutre in Meta_Data. Add
% calibratation to the appropriate structure.
if isclassfield(Meta_Data,'AFE') && ~isclassfield(Meta_Data,'epsi')
    field_name = 'AFE';
elseif isclassfield(Meta_Data,'epsi') && ~isclassfield(Meta_Data,'AFE')
    field_name = 'epsi';
elseif isclassfield(Meta_Data,'epsi') && isclassfield(Meta_Data,'AFE')
    field_name = 'epsi';
end

if isnan(str2double(Meta_Data.AFE.t1.SN)) || str2double(Meta_Data.AFE.t1.SN)==0
    fprintf('**** t1 SN (currently SN = %s)',Meta_Data.(field_name).t1.SN)
    Meta_Data.(field_name).t1.SN = input(': ','s');
end

if isnan(str2double(Meta_Data.AFE.t2.SN)) || str2double(Meta_Data.AFE.t2.SN)==0
    fprintf('**** t2 SN (currently SN = %s)',Meta_Data.(field_name).t2.SN)
    Meta_Data.(field_name).t2.SN = input(': ','s');
end

fprintf('**** t1 SN = %s ****\n',Meta_Data.(field_name).t1.SN)
fprintf('**** t2 SN = %s ****\n',Meta_Data.(field_name).t2.SN)

save(fullfile(Meta_Data.paths.data,'Meta_Data.mat'),'Meta_Data');
