function Meta_Epsi=mod_som_get_shear_probe_calibration(Meta_Epsi)

path2file1=[Meta_Epsi.shearcal_path '/Calibration_' Meta_Epsi.s1.SN '.txt'];
path2file2=[Meta_Epsi.shearcal_path '/Calibration_' Meta_Epsi.s2.SN '.txt'];

fid1=fopen(path2file1,'r');
Cal1=textscan(fid1,'%s %f %f','Delimiter',',','headerline',1);
Meta_Epsi.s1.Sv=Cal1{2}(end);

fid2=fopen(path2file2,'r');
Cal2=textscan(fid2,'%s %f %f','Delimiter',',','headerline',1);
Meta_Epsi.s2.Sv=Cal2{2}(end);

fclose(fid1);
fclose(fid2);

end