function Meta_Data=mod_som_get_shear_probe_calibration_v2(Meta_Data)

% NC edited mod_som_get_shear_probe_calibration.m May 2021
%
% Changes to Meta_Data format require a change to this function:
%   - Meta_Data no longer contains shearcal_path and shear serial numbers are
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

% shearcal_path = strrep([Meta_Data.paths.process_library,'/CALIBRATION/SHEAR_PROBES'],'//','/');
localpath=fullfile(Meta_Data.paths.process_library,'CALIBRATION','SHEAR_PROBES');

shearcal_path = strrep(localpath,'//','/');

% path2file1 = sprintf([shearcal_path '/%s/Calibration_%s.txt'], Meta_Data.(field_name).s1.SN, Meta_Data.(field_name).s1.SN);
% path2file2 = sprintf([shearcal_path '/%s/Calibration_%s.txt'], Meta_Data.(field_name).s2.SN, Meta_Data.(field_name).s2.SN);
try
    path2file1 = fullfile(shearcal_path,Meta_Data.(field_name).s1.SN,sprintf('Calibration_%s.txt',Meta_Data.(field_name).s1.SN));
    path2file2 = fullfile(shearcal_path,Meta_Data.(field_name).s2.SN,sprintf('Calibration_%s.txt',Meta_Data.(field_name).s2.SN));
catch
    path2file1 = fullfile(shearcal_path,Meta_Data.(field_name).s1.SN.',sprintf('Calibration_%s.txt',Meta_Data.(field_name).s1.SN.'));
    path2file2 = fullfile(shearcal_path,Meta_Data.(field_name).s2.SN.',sprintf('Calibration_%s.txt',Meta_Data.(field_name).s2.SN.'));
end
    

if ~strcmp(Meta_Data.(field_name).s1.SN,'000')
    try
        fid1=fopen(path2file1,'r');
        Cal1=textscan(fid1,'%s %f %f','Delimiter',',','headerline',1);
        Meta_Data.(field_name).s1.cal=Cal1{2}(end);
        fclose(fid1);
    catch err
        if strcmp(err.identifier,'MATLAB:FileIO:InvalidFid')
            warning(['Cannot find ' path2file1])
        else
            warning(['Loading ' path2file1 ' failed'])
        end
    end
end

if ~strcmp(Meta_Data.(field_name).s2.SN,'000')
    try
        fid2=fopen(path2file2,'r');
        Cal2=textscan(fid2,'%s %f %f','Delimiter',',','headerline',1);
        Meta_Data.(field_name).s2.cal=Cal2{2}(end);
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
