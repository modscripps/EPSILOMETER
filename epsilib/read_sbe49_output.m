% sbe.cal are the calibration coeficients issued after the "dcal" is sent to the SBE49
% 
% rec_ctd is a CTD records in format 0 of the SBE49 

ctd.T_raw(n_rec)  = hex2dec(rec_ctd(:,1:6));
ctd.C_raw(n_rec)  = hex2dec(rec_ctd(:,(1:6)+6));
ctd.P_raw(n_rec)  = hex2dec(rec_ctd(:,(1:6)+12));
ctd.PT_raw(n_rec) = hex2dec(rec_ctd(:,(1:4)+18));

function ctd = sbe49_ascii_get_temperature(ctd,sbe)

a0 = sbe.cal.ta0;
a1 = sbe.cal.ta1;
a2 = sbe.cal.ta2;
a3 = sbe.cal.ta3;

mv = (ctd.T_raw-524288)/1.6e7;
r = (mv*2.295e10 + 9.216e8)./(6.144e4-mv*5.3e5);
ctd.T = a0+a1*log(r)+a2*log(r).^2+a3*log(r).^3;
ctd.T = 1./ctd.T - 273.15;
return;
end

%  reads and apply calibration to the conductivity data
function ctd = sbe49_ascii_get_conductivity(ctd,sbe)
try
    g = sbe.cal.g;
    h = sbe.cal.h;
    i = sbe.cal.i;
    j = sbe.cal.j;
    tcor = sbe.cal.tcor;
    pcor = sbe.cal.pcor;
catch
    g = sbe.cal.cg;
    h = sbe.cal.ch;
    i = sbe.cal.ci;
    j = sbe.cal.cj;
    tcor = sbe.cal.ctcor;
    pcor = sbe.cal.cpcor;
end

f = ctd.C_raw/256/1000;

ctd.C = (g+h*f.^2+i*f.^3+j*f.^4)./(1+tcor.*ctd.T+pcor.*ctd.P);

return;
end

%  reads and apply calibration to the pressure data
function ctd = sbe49_ascii_get_pressure(ctd,sbe)
% ALB 04112019 Changed ctd.cal.SBEcal. to ctd.cal.
pa0 = sbe.cal.pa0;
pa1 = sbe.cal.pa1;
pa2 = sbe.cal.pa2;
ptempa0 = sbe.cal.ptempa0;
ptempa1 = sbe.cal.ptempa1;
ptempa2 = sbe.cal.ptempa2;
ptca0 = sbe.cal.ptca0;
ptca1 = sbe.cal.ptca1;
ptca2 = sbe.cal.ptca2;
ptcb0 = sbe.cal.ptcb0;
ptcb1 = sbe.cal.ptcb1;
ptcb2 = sbe.cal.ptcb2;


y = ctd.PT_raw/13107;

t = ptempa0+ptempa1*y+ptempa2*y.^2;
x = ctd.P_raw-ptca0-ptca1*t-ptca2*t.^2;
n = x*ptcb0./(ptcb0+ptcb1*t+ptcb2*t.^2);

ctd.P = (pa0+pa1*n+pa2*n.^2-14.7)*0.689476;

return;
end
