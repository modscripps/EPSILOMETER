function Meta_Epsi=get_shear_calibration(Meta_Epsi)

if (~strcmp(Meta_Epsi.s1.SN,'000') && ~strcmp(Meta_Epsi.s2.SN,'000'))
    path2file1=fullfile(Meta_Epsi.shearcal_path,Meta_Epsi.s1.SN,['Calibration_' Meta_Epsi.s1.SN '.txt']);
    path2file2=fullfile(Meta_Epsi.shearcal_path,Meta_Epsi.s2.SN,['Calibration_' Meta_Epsi.s2.SN '.txt']);
    
    fid1=fopen(path2file1,'r');
    Cal1=textscan(fid1,'%s %f %f','Delimiter',',','headerline',1);
    Meta_Epsi.s1.Sv=Cal1{2}(end);
    fclose(fid1);
    
    fid2=fopen(path2file2,'r');
    Cal2=textscan(fid2,'%s %f %f','Delimiter',',','headerline',1);
    Meta_Epsi.s2.Sv=Cal2{2}(end);    
    fclose(fid2);
end

end