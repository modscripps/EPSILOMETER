function [Meta_Data] = set_SN_temp(Meta_Data)
% Set temp probe numbers

% NC 10/7/21 - Check for 'AFE' or 'epsi' strucutre in Meta_Data. Add
% calibratation to the appropriate structure.
if isfield(obj.Meta_Data,'AFE') && ~isfield(obj.Meta_Data,'epsi')
    field_name = 'AFE';
elseif isfield(obj.Meta_Data,'epsi') && ~isfield(obj.Meta_Data,'AFE')
    field_name = 'epsi';
elseif isfield(obj.Meta_Data,'epsi') && isfield(obj.Meta_Data,'AFE')
    field_name = 'epsi';
end

fprintf('**** t1 SN (currently SN = %s)',Meta_Data.(field_name).t1.SN)
Meta_Data.(field_name).t1.SN = input(': ','s');
fprintf('**** t2 SN (currently SN = %s)',Meta_Data.(field_name).t2.SN)
Meta_Data.(field_name).t2.SN = input(': ','s');

fprintf('**** t1 SN = %s ****\n',Meta_Data.(field_name).t1.SN)
fprintf('**** t2 SN = %s ****\n',Meta_Data.(field_name).t2.SN)

save(fullfile(Meta_Data.paths.data,'Meta_Data.mat'),'Meta_Data');
