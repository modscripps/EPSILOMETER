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

% (3) add n2

% (4) add Thorpe scale

% (5) add epsilon from thorpe scale

% (6) compute gamma_chi=chi/epsilon

% (7) add ozmidov scale (epsilon/N^3)^(1/2)

% (7) add Re_b= epsilon/nu / N^2
