
function data = epsiSetup_make_empty_structure(varargin)

if nargin<1
    nSec=60*5;
elseif nargin==1
    nSec = varargin{1};
end

emptyEpsi = nan(nSec*320,1);
emptyCtd = nan(nSec*16,1);
emptyAlt = nan(nSec,1);

data.epsi.time_s = emptyEpsi;
data.epsi.dnum = emptyEpsi;
data.epsi.t1_count = emptyEpsi;
data.epsi.t2_count  = emptyEpsi;
data.epsi.s1_count = emptyEpsi;
data.epsi.s2_count = emptyEpsi;
data.epsi.a1_count = emptyEpsi;
data.epsi.a2_count = emptyEpsi;
data.epsi.a3_count = emptyEpsi;
data.epsi.t1_volt = emptyEpsi;
data.epsi.t2_volt = emptyEpsi;
data.epsi.s1_volt = emptyEpsi;
data.epsi.s2_volt = emptyEpsi;
data.epsi.a1_g = emptyEpsi;
data.epsi.a2_g = emptyEpsi;
data.epsi.a3_g = emptyEpsi;

data.ctd.time_s = emptyCtd;
data.ctd.dnum = emptyCtd;
data.ctd.P_raw = emptyCtd;
data.ctd.T_raw = emptyCtd;
data.ctd.S_raw = emptyCtd;
data.ctd.C_raw = emptyCtd;
data.ctd.PT_raw = emptyCtd;
data.ctd.P = emptyCtd;
data.ctd.T = emptyCtd;
data.ctd.S = emptyCtd;
data.ctd.C = emptyCtd;
data.ctd.sig = emptyCtd;
data.ctd.dPdt = emptyCtd;

data.alt.time_s = emptyAlt;
data.alt.dnum = emptyAlt;
data.alt.dst = emptyAlt;
