function fc_index=FPO7_cutoff(f,spec,FPO7noise)
% thcut1_t1_1_mmp
%   Usage: fc_index=thcut1_t1_1_mmp(f,spec)
%      f is a vector of frequencies
%      spec is a vector of th spectra, in volts^2/Hz
%      displ_chi_spec is a string specifying whether to plot spec
%      pr_chi is the center pressure of this spectrum
%      j in the index number of the chi estimate
%      speed is the vehicle speed, m/s
%      fc_index is the integer index of the kc, the maximum
%      frequency to which the spectrum is integrated
%   Function: determine the cut-off frequency for th spectra.
%      It is taken as the highest frequency for which a 5th-order
%      median filter exceeds the noise spectrum.  If displ_chi_spec='yes' the
%      volts^2 spectra are displayed, overlaid with the cutoff frequency and
%      the median filter.
%% Added september 14 2017 in epsilometer processing by A. LeBoyer
%  from original MMP processing 

f=f(:); spec=spec(:);
%ALB add f>0 because sometinme the frequency array have a 0; 
spec=spec(f>0);
f=f(f>0);
% fit coefficients for log10(noise spectrum) vs log10(f)
%n0=-7.8; n1=-0.0634538; n2=0.3421899; n3=-0.3141283;
%TODO change the name of the FPO7noise field so they match here
n0=FPO7noise.n0; n1=FPO7noise.n1; n2=FPO7noise.n2; n3=FPO7noise.n3;
%n0 =-0.8490;n1 = 1.5558;n2 = 0.0616;n3=-11.5939;


logf=log10(f);
noise=n0+n1.*logf+n2.*logf.^2+n3.*logf.^3;

% %TODO: why the log10(1.5)?
% logspec=log10(spec);
% medspec=medfilt1(logspec,5); % 5th order median filter
% noisy=find(medspec<log10(1.5)+noise);

%MHA
%medspec=medfilt1(spec,5); % 5th order median filter
%ALB
medspec=smoothdata(spec,'movmean',15);
% ALB really want to adjust the noise level to the high freqeuncy part (freq> 110Hz id FS is 320Hz) of
% the in-situ spectrum. The 110 Hz limit is arbitrary and means that 
% I am overlapping the noise and data spectrum between 110 and 160 Hz (if Nyquist is 160Hz).
% indhighfreq = f>110;
% ALB 09/21/2022: changing to 0.7.*f(end) becasue APEX-EPSI is sampling at
% 160 Hz and thus nyquist is 80 Hz. 0.7 is arbitrary. 
indhighfreq=find(f>0.7.*f(end));
adjust_spec=median(medspec(indhighfreq)./10.^(noise(indhighfreq)),'omitmissing');
% adjust_spec=nanmedian(medspec(indhighfreq)./10.^(noise(indhighfreq)));
% adjust_spec=1;

if adjust_spec>10
    %warning('temp noise way higher than bench')
end
SN_min=3;
%ALB I am *not* comparing the first 2 Fourier coef.
noisy=find(medspec(3:end)./adjust_spec<SN_min.*10.^(noise(3:end)));
if isempty(noisy)
	fc_index=length(f);
else
	fc_index=noisy(1);
end
if fc_index>1
	fc_index=fc_index-1;
end

