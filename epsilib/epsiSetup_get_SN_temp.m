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


possible_sensor_name={'t1','t2'};

for s=1:length(possible_sensor_name)
    wh_sensor=possible_sensor_name{s};
    if isfield(Meta_Data.(field_name),wh_sensor)

        if isnan(str2double(Meta_Data.(field_name).(wh_sensor).SN)) || str2double(Meta_Data.(field_name).(wh_sensor).SN)==0
            fprintf('**** %s SN (currently SN = %s)',wh_sensor,Meta_Data.(field_name).t1.SN)
            Meta_Data.(field_name).(wh_sensor).SN = input(': ','s');
        end
    else
        Meta_Data.(field_name).(wh_sensor).SN=nan;
        Meta_Data.(field_name).(wh_sensor).cal=nan;

    end
end

% if isnan(str2double(Meta_Data.AFE.t2.SN)) || str2double(Meta_Data.AFE.t2.SN)==0
%     fprintf('**** t2 SN (currently SN = %s)',Meta_Data.(field_name).t2.SN)
%     Meta_Data.(field_name).t2.SN = input(': ','s');
% end

fprintf('**** t1 SN = %s ****\n',Meta_Data.(field_name).t1.SN)
fprintf('**** t2 SN = %s ****\n',Meta_Data.(field_name).t2.SN)

save(fullfile(Meta_Data.paths.data,'Meta_Data.mat'),'Meta_Data');
