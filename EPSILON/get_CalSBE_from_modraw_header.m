function SBEcal=get_CalSBE_from_modraw_header(str_SBEcalcoef_header)

lines_SBEcalcoef_header=strsplit(str_SBEcalcoef_header,'\n');

%  reads and apply calibration to the conductivity data
for l=1:length(lines_SBEcalcoef_header)

    splt_linesSBE=strsplit(lines_SBEcalcoef_header{l},'=');

    switch strtrim(splt_linesSBE{1})
        case 'SERIALNO'
            SBEcal.SN=splt_linesSBE{2};
        case 'TCALDATE'
            SBEcal.TempCal_date=splt_linesSBE{2};
        case 'TA0'
            SBEcal.ta0=str2double(splt_linesSBE{2});
        case 'TA1'
            SBEcal.ta1=str2double(splt_linesSBE{2});
        case 'TA2'
            SBEcal.ta2=str2double(splt_linesSBE{2});
        case 'TA3'
            SBEcal.ta3=str2double(splt_linesSBE{2});
        case 'TOFFSET'
            SBEcal.toffset=str2double(splt_linesSBE{2});
        case 'CCALDATE'
            SBEcal.CondCal_date=splt_linesSBE{2};
        case 'CG'
            SBEcal.g=str2double(splt_linesSBE{2});
        case 'CH'
            SBEcal.h=str2double(splt_linesSBE{2});
        case 'CI'
            SBEcal.i=str2double(splt_linesSBE{2});
        case 'CJ'
            SBEcal.j=str2double(splt_linesSBE{2});
        case 'CTCOR'
            SBEcal.tcor=str2double(splt_linesSBE{2});
        case 'CPCOR'
            SBEcal.pcor=str2double(splt_linesSBE{2});
        case 'CSLOPE'
            SBEcal.pcor=str2double(splt_linesSBE{2});
        case 'PCALDATE'
            SBEcal.PresCal_date=splt_linesSBE{2};
        case 'PA0'
            SBEcal.pa0=str2double(splt_linesSBE{2});
        case 'PA1'
            SBEcal.pa1=str2double(splt_linesSBE{2});
        case 'PA2'
            SBEcal.pa2=str2double(splt_linesSBE{2});
        case 'PTCA0'
            SBEcal.ptca0=str2double(splt_linesSBE{2});
        case 'PTCA1'
            SBEcal.ptca1=str2double(splt_linesSBE{2});
        case 'PTCA2'
            SBEcal.ptca2=str2double(splt_linesSBE{2});
        case 'PTCB0'
            SBEcal.ptcb0=str2double(splt_linesSBE{2});
        case 'PTCB1'
            SBEcal.ptcb1=str2double(splt_linesSBE{2});
        case 'PTCB2'
            SBEcal.ptcb2=str2double(splt_linesSBE{2});
        case 'PTEMPA0'
            SBEcal.ptempa0=str2double(splt_linesSBE{2});
        case 'PTEMPA1'
            SBEcal.ptempa1=str2double(splt_linesSBE{2});
        case 'PTEMPA2'
            SBEcal.ptempa2=str2double(splt_linesSBE{2});
        case 'POFFSET'
            SBEcal.poffset=str2double(splt_linesSBE{2});
        otherwise
    end

end
