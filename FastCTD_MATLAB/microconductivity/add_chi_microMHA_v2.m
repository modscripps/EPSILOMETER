function FCTD=add_chi_microMHA_v2(FCTD,chi_param)
%function FCTD=add_chi_microMHA_v2(FCTD)
%Matthew Alford
%March 2019
%Version 0.2
%
%Major changes since 0.1: 
%1) chi1 and chi2 are now switched; ie chi1 is the integrated quantity and
%chi2 is stupid (time domain).
%2) Chi calc and plotting is functionized.
%3) Default chi_params structure is created and passed.
%
%Unpack microconductivity data, remove preemphasis filter, loop through in
%time windows of specified length and compute chi in those windows via a
%time-domain method and a spectral method.
%See Examine2.mlx for the core code that went into this.
%plotit=1 can be set for plots of each spectrum etc during the calculation.
%
%Dependencies:
%CenteredConv
%remove_sbe_preemphasisMHA.m
%DisplayProgress
%DoOneChi_MHA

%
%if plotting enabled:
%MySubplot
%batchelor



%%%%%%%%%%%%%%%%%%%%
%Set parameters
%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    chi_param=FCTD_DefaultChiParam;
end

chi_meta.units='degC^2/s';
chi_meta.description='chi1=stupid time-domain 6D<dTdz^2>.  chi2 is integral of wavenumber spec from kmin to kmax.';

%Create a higher-freq time vector
FCTD.microtime=linspace(FCTD.time(1),FCTD.time(end),chi_param.fs./chi_param.fslow*length(FCTD.time))';

i1=1:length(FCTD.microtime);
%Also compute fall rate in m/s
FCTD.dPdt=CenteredConv(diffs(FCTD.pressure)*chi_param.fslow,1,4*chi_param.fslow); %smooth over a few sec

%Now unpack micro and express as raw volts - convert from counts
if chi_param.fs==160
    FCTD.ucon=reshape(FCTD.uConductivity.',10*length(FCTD.time),1)/2^16;
elseif chi_param.fs==320
    FCTD.ucon=reshape(FCTD.uConductivity.',20*length(FCTD.time),1)/2^24;
end
    %reshape takes them columnwise so transpose

data=FCTD.ucon(i1); %This is the uncalibrated voltage

%This is the function that does all the work.
datac=remove_sbe_preemphasisMHA(NANinterp(data),chi_param.fs)*chi_param.gain+chi_param.offset;

%Store corrected output field with the scaling and offset applied
FCTD.ucon_corr=nan*FCTD.ucon;
FCTD.ucon_corr(i1)=datac;

%Now loop through and compute spectra for each block

%chi_param.plotit=0;
tout=FCTD.time(1):(chi_param.dt_sec/24/3600)/8:FCTD.time(end);
chi_all=nan*tout; %This is chi by simple time-domain RMS'ing, in bins
chi2_all=nan*tout; %This is chi by integrating the wavenumber spectrum
w_all=nan*tout;
disp(sprintf('Calculating chi for %7.0f scans',length(tout)))
for c=1:length(tout)

    %DisplayProgress(c,50)
    tlim=tout(c)+[-1 1].*chi_param.dt_sec/2 /24/3600; %Make a time window surrounding that center time
    
    out=DoOneChi_MHA(FCTD,tlim,chi_param);
    
    w_all(c)=out.w;
    chi2_all(c)=out.chi_stupid;
    chi_all(c)=out.chi1;
    
    % Print progress status
    if mod(c,20)==0 && mod(c,100)~=0
        fprintf(num2str(c))
    elseif mod(c,100)==0
        fprintf([num2str(c) '\n'])
    elseif c==length(tout)
        fprintf('. \n')
    else
        fprintf('.')
    end
    
end



%% Now output to the FCTD grid
FCTD.chi=interp1(tout,chi_all,FCTD.time);
FCTD.chi2=interp1(tout,chi2_all,FCTD.time);
FCTD.w=interp1(tout,w_all,FCTD.time);

FCTD.chi_param=chi_param;
FCTD.chi_meta=chi_meta;



