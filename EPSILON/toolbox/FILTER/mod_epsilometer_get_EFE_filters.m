function H=mod_epsilometer_get_EFE_filters(Meta_Data,f)
%function H=get_filters_MADRE(Meta_Data,f)
%% Define the electronic filter, ADC filters.
%  The electronics filter Helec = charge amp filter (Hca) + Gain (Hg)
%  The ADC filter is a simple sinc^4 and the gain of the ADC is 2. Other options are available. 
%
%  TODO: enter argument to be able to change the ADC filter without hard
%  coding
% Feb 2019 ALB

switch Meta_Data.epsi.s1.ADCfilter
    case 'sinc4'
        Hs1filter=(sinc(f/(2*f(end)))).^4;
end
switch Meta_Data.epsi.t1.ADCfilter
    case 'sinc4'
        Ht1filter=(sinc(f/(2*f(end)))).^4;
end
switch Meta_Data.epsi.a1.ADCfilter
    case 'sinc4'
        Ha1filter=(sinc(f/(2*f(end)))).^4;
end

% shear channels
%charge amp filter

%ca_filter = load('FILTER/charge_coeffilt.mat');
% calibrator gain 2.5 dB ~ Vout/Vin=1.33 for a shear probe of 750pF

%ca_filter = load('FILTER/charge_coeffilt_09312019.mat'); %from network analysis
% ca_filter = load('FILTER/charge_coeffilt.mat'); %from network analysis

% ALB until now I never used the AA TF. May 30 th 2020
ca_filter = load('FILTER/shr_chramp_waa.mat'); %from spice WITH AA filter. 
epsi_ca   = interp1(ca_filter.freq,ca_filter.coef_filt ,f);
gain_ca   = 1; %try with a 0 dB gain.
H.electshear= epsi_ca*gain_ca;% 
H.gainshear=1;
H.adcshear=H.gainshear.* Hs1filter;
H.shear=(H.electshear .* H.adcshear).^2;

%% FPO7 channels

H.gainFPO7=1;
H.electFPO7 = H.gainFPO7.*Ht1filter;
%speed% convert to m/s
%tau=0.005 * speed^(-0.32); % thermistor time constant
H.magsq=@(speed)(1 ./ (1+((2*pi*(0.005 * speed^(-0.32))).*f).^2)); % magnitude-squared no units
H.phase=@(speed)(-2*atan( 2*pi*f*(0.005 * speed^(-0.32))));   % no units
switch Meta_Data.Hardware.EFE.temperature
    case 'Tdiff'
        Tdiff_filter = load('FILTER/Tdiff_filt.mat');
        Tdiff_H = interp1(Tdiff_filter.freq,Tdiff_filter.coef_filt ,f);
        H.Tdiff=Tdiff_H;
        H.FPO7=@(speed)(H.electFPO7.^2 .* H.magsq(speed) .* H.Tdiff.^2);
    otherwise
        H.FPO7=@(speed)(H.electFPO7.^2 .* H.magsq(speed));
end



%% Accel channels
H.gainAccel  = 1;
H.electAccel = (H.gainAccel.*Ha1filter).^2;

end

