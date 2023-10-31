function chi_param=FCTD_DefaultChiParam 

%Define sample freqs
fslow=16;
fs=160;

%Offset/gain - iteratively plot below and adjust values here.
% These are the gain and offsets determined from matching the low-freq
% cond.
%cond_corr=gain*voltage + offset;

gain=22;%1/2900;
offset=2.72;%not important for chi calc
dTdC=10;% in degC/(S/m).  Nominal value.

nsec=0.5; %number of secs going into the FFT window.  First tried with 4.
dt_sec=2; %length of the data window.  F, dirst tried with 20.

kmin=2;kmax=200; %min, max wavenumbers in cpm for integration
min_spd=0.5; %min speed to try calculation

D=1.4e-7; %thermal diffusivity

chi_param.fslow=fslow;
chi_param.fs=fs;
chi_param.gain=gain;
chi_param.offset=offset;
chi_param.dTdC=dTdC;
chi_param.nsec=nsec;
chi_param.dt_sec=2;
chi_param.kmin=2;
chi_param.kmax=200;
chi_param.min_spd=0.5;
chi_param.D=D;
chi_param.mfile='add_chi_microMHA_v2.m';
chi_param.version='0.2';
chi_param.plotit=0;
