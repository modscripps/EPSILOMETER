function SBEcal=get_CalSBE_from_modraw(str)
% %  reads and apply calibration to the conductivity data
% 
%     fid=fopen(filename);
%     line=fgetl(fid);
%     SBEcal.SN=line(end-4:end);
%     % Temperature Cal
%     line=fgetl(fid);
%     SBEcal.TempCal_date=line(end-10:end);
%     line=fgetl(fid);
%     SBEcal.ta0=str2double(line(end-12:end));
%     line=fgetl(fid);
%     SBEcal.ta1=str2double(line(end-12:end));
%     line=fgetl(fid);
%     SBEcal.ta2=str2double(line(end-13:end));
%     line=fgetl(fid);
%     SBEcal.ta3=str2double(line(end-12:end));
%     line=fgetl(fid);
%     SBEcal.toffset=str2double(line(end-12:end));
%     % Conductivity Cal
%     line=fgetl(fid);
%     SBEcal.CondCal_date=line(end-10:end);
%     line=fgetl(fid);
%     SBEcal.g=str2double(line(end-12:end));
%     line=fgetl(fid);
%     SBEcal.h=str2double(line(end-12:end));
%     line=fgetl(fid);
%     SBEcal.i=str2double(line(end-12:end));
%     line=fgetl(fid);
%     SBEcal.j=str2double(line(end-12:end));
%     line=fgetl(fid);
%     SBEcal.pcor=str2double(line(end-12:end));
%     line=fgetl(fid);
%     SBEcal.tcor=str2double(line(end-12:end));
%     line=fgetl(fid);
%     SBEcal.cslope=str2double(line(end-12:end));
%     % Pressure Cal
%     line=fgetl(fid);
%     SBEcal.PresCal_date=line(end-10:end);
%     line=fgetl(fid);
%     SBEcal.pa0=str2double(line(end-12:end));
%     line=fgetl(fid);
%     SBEcal.pa1=str2double(line(end-12:end));
%     line=fgetl(fid);
%     SBEcal.pa2=str2double(line(end-12:end));
%     line=fgetl(fid);
%     SBEcal.ptca0=str2double(line(end-12:end));
%     line=fgetl(fid);
%     SBEcal.ptca1=str2double(line(end-12:end));
%     line=fgetl(fid);
%     SBEcal.ptca2=str2double(line(end-12:end));
%     line=fgetl(fid);
%     SBEcal.ptcb0=str2double(line(end-12:end));
%     line=fgetl(fid);
%     SBEcal.ptcb1=str2double(line(end-12:end));
%     line=fgetl(fid);
%     SBEcal.ptcb2=str2double(line(end-12:end));
%     line=fgetl(fid);
%     SBEcal.ptempa0=str2double(line(end-12:end));
%     line=fgetl(fid);
%     SBEcal.ptempa1=str2double(line(end-12:end));
%     line=fgetl(fid);
%     SBEcal.ptempa2=str2double(line(end-12:end));
%     line=fgetl(fid);
%     SBEcal.poffset=str2double(line(end-12:end));
%     fclose(fid);


%  reads and apply calibration to the conductivity data

    fid=fopen(filename);
    line=fgetl(fid);
    SBEcal.SN=line(end-4:end);
    % Temperature Cal
    line=fgetl(fid);
    SBEcal.TempCal_date=line(end-10:end);
    line=fgetl(fid);
    SBEcal.ta0=str2double(line(strfind(line,'=')+1:end));
    line=fgetl(fid);
    SBEcal.ta1=str2double(line(strfind(line,'=')+1:end));
    line=fgetl(fid);
    SBEcal.ta2=str2double(line(strfind(line,'=')+1:end));
    line=fgetl(fid);
    SBEcal.ta3=str2double(line(strfind(line,'=')+1:end));
    line=fgetl(fid);
    SBEcal.toffset=str2double(line(strfind(line,'=')+1:end));
    % Conductivity Cal
    line=fgetl(fid);
    SBEcal.CondCal_date=line(end-10:end);
    line=fgetl(fid);
    SBEcal.g=str2double(line(strfind(line,'=')+1:end));
    line=fgetl(fid);
    SBEcal.h=str2double(line(strfind(line,'=')+1:end));
    line=fgetl(fid);
    SBEcal.i=str2double(line(strfind(line,'=')+1:end));
    line=fgetl(fid);
    SBEcal.j=str2double(line(strfind(line,'=')+1:end));
    line=fgetl(fid);
    SBEcal.pcor=str2double(line(strfind(line,'=')+1:end));
    line=fgetl(fid);
    SBEcal.tcor=str2double(line(strfind(line,'=')+1:end));
    line=fgetl(fid);
    SBEcal.cslope=str2double(line(strfind(line,'=')+1:end));
    % Pressure Cal
    line=fgetl(fid);
    SBEcal.PresCal_date=line(end-10:end);
    line=fgetl(fid);
    SBEcal.pa0=str2double(line(strfind(line,'=')+1:end));
    line=fgetl(fid);
    SBEcal.pa1=str2double(line(strfind(line,'=')+1:end));
    line=fgetl(fid);
    SBEcal.pa2=str2double(line(strfind(line,'=')+1:end));
    line=fgetl(fid);
    SBEcal.ptca0=str2double(line(strfind(line,'=')+1:end));
    line=fgetl(fid);
    SBEcal.ptca1=str2double(line(strfind(line,'=')+1:end));
    line=fgetl(fid);
    SBEcal.ptca2=str2double(line(strfind(line,'=')+1:end));
    line=fgetl(fid);
    SBEcal.ptcb0=str2double(line(strfind(line,'=')+1:end));
    line=fgetl(fid);
    SBEcal.ptcb1=str2double(line(strfind(line,'=')+1:end));
    line=fgetl(fid);
    SBEcal.ptcb2=str2double(line(strfind(line,'=')+1:end));
    line=fgetl(fid);
    SBEcal.ptempa0=str2double(line(strfind(line,'=')+1:end));
    line=fgetl(fid);
    SBEcal.ptempa1=str2double(line(strfind(line,'=')+1:end));
    line=fgetl(fid);
    SBEcal.ptempa2=str2double(line(strfind(line,'=')+1:end));
    line=fgetl(fid);
    SBEcal.poffset=str2double(line(strfind(line,'=')+1:end));
    fclose(fid);

end
