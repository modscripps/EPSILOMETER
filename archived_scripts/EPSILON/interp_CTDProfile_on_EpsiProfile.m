function Profile=interp_CTDProfile_on_EpsiProfile(CTDProfile,EpsiProfile)

rem_nan=@(x) (fillmissing(x,'linear'));

Profile=EpsiProfile;
Profile.P=rem_nan(interp1(rem_nan(CTDProfile.ctdtime),CTDProfile.P,EpsiProfile.epsitime));
Profile.P=filloutliers(Profile.P,'center','movmedian',1000);
Profile.T=rem_nan(interp1(rem_nan(CTDProfile.ctdtime),CTDProfile.T,EpsiProfile.epsitime));
Profile.S=rem_nan(interp1(rem_nan(CTDProfile.ctdtime),CTDProfile.S,EpsiProfile.epsitime));
Profile.time=EpsiProfile.epsitime;

