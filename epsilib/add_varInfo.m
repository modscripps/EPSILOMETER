function [Profile] = add_varInfo(Profile);

% The order of the fields in varInfo will define the order of fields in
% Profile, so make sure there are the same number of fields!
if isfield(Profile.ctd,'dnum')
    Profile.varInfo.ctd.dnum = {'datenum','Matlab datenum'};
    Profile.varInfo.ctd.time_s = {'time','seconds since Jan 1 1970'};
else
    Profile.varInfo.ctd.time_s = {'time','seconds since power on'};
end
Profile.varInfo.ctd.P = {'CTD P','db'};
Profile.varInfo.ctd.dPdt = {'CTD diff(P)/diff(ctdtime)','db s^{-1}'};
Profile.varInfo.ctd.T = {'CTD temperature','^\circC'};
Profile.varInfo.ctd.C = {'CTD conductivity','S/m'};
Profile.varInfo.ctd.S = {'CTD salinity','psu'};
Profile.varInfo.th = {'Potential temperature','^\circC'};
Profile.varInfo.ctd.sgth = {'CTD potential density (sigma-theta)',''};

Profile.varInfo.epsi.a1_g = {'acceleration sensor 1 timeseries in this scan','[g]'};
Profile.varInfo.epsi.a2_g = {'acceleration sensor 2 timeseries in this scan','[g]'};
Profile.varInfo.epsi.a3_g = {'acceleration sensor 3 timeseries in this scan','[g]'};
Profile.varInfo.epsi.s1_volt = {'shear sensor 1 timeseries in this scan','Volts'};
Profile.varInfo.epsi.s2_volt = {'shear sensor 2 timeseries in this scan','Volts'};
Profile.varInfo.epsi.t1_volt = {'temperature sensor 1 timeseries in this scan','Volts'};
Profile.varInfo.epsi.t2_volt = {'temperature sensor 2 timeseries in this scan','Volts'};
if isfield(Profile.epsi,'c_count')
Profile.varInfo.epsi.c_count = {'additional sensor timseries in the scan' 'count'};
end

Profile.varInfo.nbscan = {'number of scans','#'};
Profile.varInfo.nfft = {'',''};
Profile.varInfo.nfftc = {'',''};
Profile.varInfo.tscan = {'length of scan window','s'};
Profile.varInfo.fpump = {'',''};
Profile.varInfo.dnum = {'datenum','Matlab datenum'};
Profile.varInfo.pr = {'CTD pressure','db'};
Profile.varInfo.w = {'fall speed','m s^{-1}'};
Profile.varInfo.t = {'temperature','^\circC'};
Profile.varInfo.s = {'salinity','psu'};
Profile.varInfo.kvis = {'kinematic viscosity',''};
Profile.varInfo.epsilon = {'turbulent kinetic energy dissipation rate calculated from Ps_shear_k', ''};
Profile.varInfo.epsilon_co = {'turbulent kinetic energy dissipation rate calculated from Ps_shear_co_k', ''};
Profile.varInfo.chi = {'temperature gradient dissipation rate','°C^2 s^{-1}'};
Profile.varInfo.sh_fc = {'shear cutoff frequency, 1=uncorrected, 2=coherence-corrected', 'Hz'};
Profile.varInfo.tg_fc = {'temperature gradient cutoff frequency, 1=uncorrected, 2=coherence-corrected', 'Hz'};
Profile.varInfo.flag_tg_fc = {'temperature gradient cut off frequency is very high','0/1'};
Profile.varInfo.ind_range_ctd = {'1st and last indices of CTD arrays in this Profile that go into each scan window',''};
Profile.varInfo.ind_range_epsi = {'1st and last indices of Epsi arrays in this Profile that go into each scan window',''};
Profile.varInfo.f = {'frequency','Hz'};
Profile.varInfo.k = {'wavenumber','cycles m^-^1'};
Profile.varInfo.Pa_g_f = {'accleration frequency power spectrum', '[g]^2 Hz^{-1}'};
Profile.varInfo.Ps_volt_f = {'shear frequency power spectrum', 'Volts^2 Hz^{-1}'};
Profile.varInfo.Ps_velocity_f = {'Velocity frequency power spectrum', 'm^2 s^{-2} Hz^{-1}'};
Profile.varInfo.Ps_shear_k = {'shear wavenumber power spectrum', 's{-1} cpm^{-1}'};
Profile.varInfo.Ps_shear_co_k = {'coherence-corrected shear frequency power spectrum (full profile coherence with a3 channel has been removed)', ''};
Profile.varInfo.Pt_volt_f = {'temperature frequency power spectrum','Volts^2 Hz{-1}'};
Profile.varInfo.Pt_Tg_k = {'temperature gradient wavenumber power spectrum', 'C^2 s{-1} cpm^{-1}'};
Profile.varInfo.Cs1a3_full = {'coherence betwen s1 and a3 channels between Meta_Data.PROCESS.Prmin and Meta_Data.PROCESS.Prmax',''};
Profile.varInfo.Cs2a3_full = {'coherence betwen s2 and a3 channels between Meta_Data.PROCESS.Prmin and Meta_Data.PROCESS.Prmax',''};

end
