% epsi_profile_class
%
% Currently just a dumping ground for functions we want to add at the
% Profile level

% (1) get_hab - get height above bottom
if isfield(Profile,'alt')
    ctdZ = interp1(Profile.ctd.time_s,Profile.ctd.z,Profile.alt.time_s);
    hab = Profile.alt.dst;
    hab(hab>35) = nan;
end

% (2) Interpolate ctd and alt data to epsi timestamp