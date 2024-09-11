function obj = epsiSetup_make_empty_structure(obj,varargin)

if nargin<2
    nSec=60*2; %2 minutes
elseif nargin==2
    nSec = varargin{1};
end

% Make some empty obj strucutures approximating their relative size. These
% don't have to be totally accurate. 
emptyEpsi = nan(nSec*320,1); % epsi samples at 320 Hz
emptyCtd  = nan(nSec*16,1); % ctd samples at 16 Hz
emptyAlt  = nan(nSec,1); % alt samples around 1 Hz
emptyVnav = nan(nSec*40,1); %vnav samples around 40 Hz
emptyGps  = nan(nSec,1); % alt samples around 1 Hz
emptyTtv  = nan(nSec*20,1); %Ti ttv samples around 20 Hz

obj.epsi.time_s = emptyEpsi;
obj.epsi.dnum = emptyEpsi;
% obj.epsi.t1_count = emptyEpsi;
% obj.epsi.t2_count  = emptyEpsi;
% obj.epsi.s1_count = emptyEpsi;
% obj.epsi.s2_count = emptyEpsi;
% obj.epsi.a1_count = emptyEpsi;
% obj.epsi.a2_count = emptyEpsi;
% obj.epsi.a3_count = emptyEpsi;
obj.epsi.t1_volt = emptyEpsi;
obj.epsi.t2_volt = emptyEpsi;
obj.epsi.s1_volt = emptyEpsi;
obj.epsi.s2_volt = emptyEpsi;
obj.epsi.a1_g = emptyEpsi;
obj.epsi.a2_g = emptyEpsi;
obj.epsi.a3_g = emptyEpsi;

obj.ctd.time_s = emptyCtd;
obj.ctd.dnum = emptyCtd;
obj.ctd.P_raw = emptyCtd;
obj.ctd.T_raw = emptyCtd;
obj.ctd.S_raw = emptyCtd;
obj.ctd.C_raw = emptyCtd;
obj.ctd.PT_raw = emptyCtd;
obj.ctd.P = emptyCtd;
obj.ctd.T = emptyCtd;
obj.ctd.S = emptyCtd;
obj.ctd.C = emptyCtd;
obj.ctd.sig = emptyCtd;
obj.ctd.dPdt = emptyCtd;
obj.ctd.z = emptyCtd;
obj.ctd.dzdt = emptyCtd;

obj.alt.time_s = emptyAlt;
obj.alt.dnum = emptyAlt;
obj.alt.dst = emptyAlt;

obj.vnav.time_s       = emptyVnav;
obj.vnav.dnum         = emptyVnav;
obj.vnav.compass      = repmat(emptyVnav,1,3);
obj.vnav.acceleration = repmat(emptyVnav,1,3);
obj.vnav.gyro         = repmat(emptyVnav,1,3);

obj.ttv.time_s    = emptyTtv;
obj.ttv.dnum      = emptyTtv;
obj.ttv.dtof     = emptyTtv;
obj.ttv.tof_up   = emptyTtv;
obj.ttv.tof_down = emptyTtv;
obj.ttv.vfr      = emptyTtv;


obj.gps.dnum      = emptyGps;
obj.gps.latitude  = emptyGps;
obj.gps.longitude = emptyGps;

