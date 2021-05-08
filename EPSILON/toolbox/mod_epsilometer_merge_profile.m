function Profile=mod_epsilometer_merge_profile(Meta_Data, CTDProfile,EpsiProfile,Prmin,Prmax)
%  Profile=mod_epsilometer_merge_profile(CTDProfile,EpsiProfile,Prmin,Prmax)

Profile=structfun(@(x) x(CTDProfile.P>=Prmin & CTDProfile.P<=Prmax),CTDProfile,'un',0);

for fields=fieldnames(EpsiProfile)'
    Profile.(fields{1})=EpsiProfile.(fields{1}) ...
         (EpsiProfile.epsitime>=CTDProfile.ctdtime(1) & ...
          EpsiProfile.epsitime<=CTDProfile.ctdtime(end));
end

% Get rid of some of the fields we don't need in L1 Profile
if isfield(Profile,'aux1time') && all(Profile.aux1time==Profile.ctdtime)
    Profile = rmfield(Profile,'aux1time');
end
if isfield(Profile,'flagSDSTR')
    Profile = rmfield(Profile,'flagSDSTR');
end
if isfield(Profile,'ramp_count')
    Profile = rmfield(Profile,'ramp_count');
end
