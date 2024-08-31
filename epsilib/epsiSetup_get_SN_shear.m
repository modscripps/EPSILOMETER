function [Meta_Data] = epsiSetup_get_SN_shear(Meta_Data)
% Set shear probe number and get calibration data

% NC 10/7/21 - Check for 'AFE' or 'epsi' strucutre in Meta_Data. Add
% calibratation to the appropriate structure.
% if isclassfield(obj.Meta_Data,'AFE') && ~isclassfield(obj.Meta_Data,'epsi')
%     field_name = 'AFE';
% elseif isclassfield(obj.Meta_Data,'epsi') && ~isclassfield(Meta_Data,'AFE')
%     field_name = 'epsi';
% elseif isclassfield(Meta_Data,'epsi') && isclassfield(Meta_Data,'AFE')
%     field_name = 'epsi';
% end

if isclassfield(Meta_Data,'AFE') && ~isclassfield(Meta_Data,'epsi')
    field_name = 'AFE';
elseif isclassfield(Meta_Data,'epsi') && ~isclassfield(Meta_Data,'AFE')
    field_name = 'epsi';
elseif isclassfield(Meta_Data,'epsi') && isclassfield(Meta_Data,'AFE')
    field_name = 'epsi';
end

possible_sensor_name={'s1','s2'};

for s=1:length(possible_sensor_name)
    wh_sensor=possible_sensor_name{s};
    if isfield(Meta_Data.(field_name),wh_sensor)
        % get SN
        if isnan(str2double(Meta_Data.(field_name).(wh_sensor).SN)) || str2double(Meta_Data.(field_name).(wh_sensor).SN)==0
            fprintf('**** %s SN (currently SN = %s, cal = %3.2f)', ...
                wh_sensor,...
                Meta_Data.(field_name).(wh_sensor).SN, ...
                Meta_Data.(field_name).(wh_sensor).cal)
            Meta_Data.(field_name).(wh_sensor).SN = input(': ','s');
        end
    else
        Meta_Data.(field_name).(wh_sensor).SN=nan;
        Meta_Data.(field_name).(wh_sensor).cal=nan;
    end
end

% if isnan(str2double(Meta_Data.(field_name).s2.SN)) || str2double(Meta_Data.(field_name).s2.SN)==0
%     fprintf('**** s2 SN (currently SN = %s, cal = %3.2f)',Meta_Data.(field_name).s2.SN,Meta_Data.(field_name).s2.cal)
%     Meta_Data.(field_name).s2.SN = input(': ','s');
% end

Meta_Data=mod_som_get_shear_probe_calibration_v2(Meta_Data);

fprintf('**** s1 SN = %s, cal = %3.2f ****\n',Meta_Data.(field_name).s1.SN,Meta_Data.(field_name).s1.cal)
fprintf('**** s2 SN = %s, cal = %3.2f ****\n',Meta_Data.(field_name).s2.SN,Meta_Data.(field_name).s2.cal)

save(fullfile(Meta_Data.paths.data,'Meta_Data.mat'),'Meta_Data');
