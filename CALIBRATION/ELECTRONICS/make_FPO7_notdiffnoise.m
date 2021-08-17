% Make FPO7_notdiffnoise.mat
% 
% Nicole Couto | June 2021
% -------------------------------------------------------------------------
% Noise floor is Johnson noise + Amplifier noise
%
% Johnson noise = vn^2/Nquist 
%
%   Nyquist = (1/2)*325 Hz
%   vn^2 = 4*kB*T*R
%
%   kB = 1.38e-23 J/K
%   R  = 150 kOhm 
%   T  = 273 + 20 K

f_epsi=.01:.01:160;
Nyf = 325/2;

kB = 1.38e-23;
R = 150e3;
T = 293;
vn2 = 4*kB*T*R;
noise_johnson = vn2/Nyf * ones(size(f_epsi)); 

% Amplifier noise = Input voltage noise + Input current noise, run through
% the Tdiff response
%
% The voltage noise at 20 Hz is 100e-9 A/sqrt(Hz).
% The current noise at 20 Hz is 3e-12 A/sqrt(Hz). Multiply that by the
% bridge resistance of 150 kOhm / 2.
% The total amplifier noise is the sum of the squared voltage and current noise

noise_voltage = 100e-9;
noise_current = 3e-12*R/2;
noise_amp = noise_current.^2 + noise_voltage.^2;

% Then this gets multiplied by 20/f (1/f that is unity at 20 Hz) and run
% through the Tdiff filter.
% Instead of the Tdiff filter, we use a fudge factor, ff=200.

ff = 200; %fudge factor
noise_amp = (noise_amp.*(20./f_epsi));

fpo7noise = noise_johnson +ff*noise_amp;

% Now, turn this noise into coefficients that will be called in later
% functions to create a polynomial fit
% noise = n0+n1.*logf+n2.*logf.^2+n3.*logf.^3;

logf = log10(f_epsi);
F = fit(logf(:),log10(fpo7noise(:)),'poly1');
% F(x) = p1*x + p2
n0 = F.p2;
n1 = F.p1;
n2 = 0;
n3 = 0;

save('/Users/ncouto/GitHub/EPSILOMETER/CALIBRATION/ELECTRONICS/FPO7_notdiffnoise_new',...
    'n0','n1','n2','n3')
