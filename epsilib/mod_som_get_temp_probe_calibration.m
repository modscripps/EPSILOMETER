function Meta_Data=mod_som_get_temp_probe_calibration(Meta_Data)

%
% Changes to Meta_Data format require a change to this function:
%   - tempcal_path and temp serial numbers are
%     now stored in Meta_Data.AFE

% NC 10/7/21 - Check for 'AFE' or 'epsi' strucutre in Meta_Data. Add
% calibratation to the appropriate structure.
if isfield(Meta_Data,'AFE') && ~isfield(Meta_Data,'epsi')
    field_name = 'AFE';
elseif isfield(Meta_Data,'epsi') && ~isfield(Meta_Data,'AFE')
    field_name = 'epsi';
elseif isfield(Meta_Data,'epsi') && isfield(Meta_Data,'AFE')
    field_name = 'epsi';
end

% tempcal_path = strrep([Meta_Data.paths.process_library,'/CALIBRATION/SHEAR_PROBES'],'//','/');
% TODO change local file to Meta_Data.AFE....
localpath=fullfile(Meta_Data.paths.process_library,'CALIBRATION','FPO7');

tempcal_path = strrep(localpath,'//','/');

try
    path2file1 = fullfile(tempcal_path,Meta_Data.(field_name).t1.SN,sprintf('Calibration_%s.txt',Meta_Data.(field_name).t1.SN));
    path2file2 = fullfile(tempcal_path,Meta_Data.(field_name).t2.SN,sprintf('Calibration_%s.txt',Meta_Data.(field_name).t2.SN));
catch
    path2file1 = fullfile(tempcal_path,Meta_Data.(field_name).t1.SN.',sprintf('Calibration_%s.txt',Meta_Data.(field_name).t1.SN.'));
    path2file2 = fullfile(tempcal_path,Meta_Data.(field_name).t2.SN.',sprintf('Calibration_%s.txt',Meta_Data.(field_name).t2.SN.'));
end
    

if ~strcmp(Meta_Data.(field_name).t1.SN,'000')
    try
        fid1=fopen(path2file1,'r');
        Cal1=textscan(fid1,'%s %f %f','Delimiter',',','headerline',1);
        Meta_Data.(field_name).t1.cal=Cal1{2}(end);
        Meta_Data.(field_name).t1.dTdV=Cal1{2}(end);
        fclose(fid1);
    catch err
        if strcmp(err.identifier,'MATLAB:FileIO:InvalidFid')
            warning(['Cannot find ' path2file1])
        else
            warning(['Loading ' path2file1 ' failed'])
        end
    end
end

if ~strcmp(Meta_Data.(field_name).t2.SN,'000')
    try
        fid2=fopen(path2file2,'r');
        Cal2=textscan(fid2,'%s %f %s','Delimiter',',','headerline',1);
        Meta_Data.(field_name).t2.cal=Cal2{2}(end);
        Meta_Data.(field_name).t2.dTdV=Cal2{2}(end);
        fclose(fid2);
    catch err
        if strcmp(err.identifier,'MATLAB:FileIO:InvalidFid')
            warning(['Cannot find ' path2file2])
        else
            warning(['Loading ' path2file2 ' failed'])
        end
    end
end
end
