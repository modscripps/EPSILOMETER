function Meta_Data=mod_som_get_shear_probe_calibration_v2(Meta_Data)

% NC edited mod_som_get_shear_probe_calibration.m May 2021
%
% Changes to Meta_Data format require a change to this function:
%   - Meta_Data no longer contains shearcal_path and shear serial numbers are
%     now stored in Meta_Data.AFE

shearcal_path = strrep([Meta_Data.processpath,'/CALIBRATION/SHEAR_PROBES'],'//','/');

path2file1 = sprintf([shearcal_path '/%s/Calibration_%s.txt'], Meta_Data.AFE.s1.SN, Meta_Data.AFE.s1.SN);
path2file2 = sprintf([shearcal_path '/%s/Calibration_%s.txt'], Meta_Data.AFE.s2.SN, Meta_Data.AFE.s2.SN);

fid1=fopen(path2file1,'r');
Cal1=textscan(fid1,'%s %f %f','Delimiter',',','headerline',1);
Meta_Data.AFE.s1.cal=Cal1{2}(end);

fid2=fopen(path2file2,'r');
Cal2=textscan(fid2,'%s %f %f','Delimiter',',','headerline',1);
Meta_Data.AFE.s2.cal=Cal2{2}(end);

fclose(fid1);
fclose(fid2);

end