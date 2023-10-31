function [coeff_dTdC,coeff_uC] = FCTD_ucon_cal(fctd,chi_param)
%FCTD_ucon_calcoeff calculates the calibration coefficents for temperature
%gradient and micro-conducitivity needed for processing microconductivity data
% -> dTdC and gain between conductivity and micro-conductivity 
%
%Bethan Wynne-Cattanach
%June 2020
%
%
%  fctd is the structure containing the full time series
%
%

%This calculates the calibration coefficient using spectra of 
%temp gradient, conductivity and micro-conductivity.
% For efficinecy it picks depths evenly spaced by 200m from 100m to the max
%depth. Contains an if statement to catch any segments that are when the
%fctd is changing direction.
if nargin < 2
    chi_param=FCTD_DefaultChiParam;
end

if max(fctd.pressure)>100
% max_d=round(max(fctd.pressure,[],'all')); %max depth rounded to closest 100m
% min_d=round(min(fctd.pressure,[],'all')); 
% depths=[min_d:5:max_d];
% 
% indx=[];
% for i=1:length(depths)
%     temp=find(abs(fctd.pressure-depths(i))<0.02);
%     indx=[indx; temp];
% end

%Look for the middle of the profiles
mins=find(islocalmin(fctd.pressure));
maxs=find(islocalmax(fctd.pressure));
if ~isempty(mins) && ~isempty(maxs)
    if length(mins)>length(maxs)
        maxs=[maxs; length(fctd.pressure)];
    elseif length(mins)<length(maxs)
        mins=[mins; length(fctd.pressure)];
    end
    indx=round((maxs+mins)/2);
else 
    max_d=round(max(fctd.pressure,[],'all')); %max depth rounded to closest 100m
    min_d=round(min(fctd.pressure,[],'all')); 
    depths=min_d:100:max_d;
    [~,indx]=closest(fctd.pressure,depths);
end

indx=sort(indx);
if chi_param.plotit
figure
end
for i=1:length(indx)
    c=indx(i);
    tlim=fctd.time(c)+[-50 50].*(1/24/3600);

    i1=find(fctd.time>tlim(1)&fctd.time<tlim(2));
    if length(i1)<160
        coeff_dTdC(i)=NaN;
        coeff_uC(i)=NaN;
    else
        %i2=find(fctd.microtime>tlim(1)&fctd.microtime<tlim(2));

        % if all(diff(fctd.pressure(i1))>0) || all(diff(fctd.pressure(i1))<0)
    %     if mod(length(i2),2)==1
    %         i2=[i2; i2(end)+1];
    %     end

        % Calculate spectra
        time_sec=(fctd.time-fctd.time(1))*3600*24;
        microtime_sec=(fctd.microtime-fctd.microtime(1))*3600*24;

        %conductivity
        [f1c,P1c] = myspectrum(NANinterp(fctd.conductivity(i1)),160,1/16,0); %Set NFFT equal to window length    
        [f2c,P2c] = myspectrum(diff(NANinterp(fctd.conductivity(i1)))*chi_param.fslow,160,1/16,0);


        %[f_c,amp_c]=spectrum_BLWC(time_sec(i1),NANinterp(f ctd.conductivity(i1)),'seglength',160,'overlap','yes','window','hamming');
        %[f_cdiff,amp_cdiff]=spectrum_BLWC(time_sec(i1(1:end-1)),diff(NANinterp(fctd.conductivity(i1))),'seglength',160,'overlap','yes','window','hamming');  
        %temperature
        [f1t,P1t] = myspectrum(NANinterp(fctd.temperature(i1)),160,1/16,0); %Set NFFT equal to window length    
        [f2t,P2t] = myspectrum(diff(NANinterp(fctd.temperature(i1)))*chi_param.fslow,160,1/16,0);

         %[f_t,amp_t]=spectrum_BLWC(time_sec(i1),NANinterp(fctd.temperature(i1)),'seglength',160,'overlap','yes','window','hamming');
         %[f_tdiff,amp_tdiff]=spectrum_BLWC(time_sec(i1(1:end-1)),diff(NANinterp(fctd.temperature(i1))),'seglength',160,'overlap','yes','window','hamming');
        %microconductivity
        %[f_uc,amp_uc]=spectrum_BLWC(microtime_sec(i2),NANinterp(remove_sbe_preemphasisMHA(fctd.ucon(i2),chi_param.fs)),'seglength',1600,'window','hamming');
       [f1uc,P1uc]=myspectrum(remove_sbe_preemphasisMHA(NANinterp(fctd.ucon(fctd.microtime>tlim(1)&fctd.microtime<tlim(2))),chi_param.fs),1600,1/160,0);
       % [f1uc,P1uc]=myspectrum(NANinterp(fctd.ucon(fctd.microtime>tlim(1)&fctd.microtime<tlim(2))),1600,1/160,0);

    %     coeff_dTdC(i)=sqrt(nanmean(amp_tdiff(1:20))./nanmean(amp_cdiff(1:20)));
    %     coeff_uC(i)=sqrt(nanmean(amp_c)./nanmean(amp_uc(1:length(amp_c))));
        coeff_dTdC(i)=sqrt(nanmean(P2t(1:20))./nanmean(P2c(1:20)));
        coeff_uC(i)=sqrt(nanmean(P1c(1:20))./nanmean(P1uc(1:20)));

    
        if chi_param.plotit
            a1=subplot(3,1,1);
            loglog(f1uc,P1uc*coeff_uC(i)^2)
            hold on
            loglog(f1c,P1c)
            loglog(f1uc,P1uc)
            ylabel('PSD (S^2m^{-2}Hz^{-1})')
            legend('\mu-cond. (post calibration)','cond.','\mu-cond. (pre calibration)')
            title(sprintf('gain=%.1f',coeff_uC(i)))
            grid on
            a2=subplot(3,1,2);
            loglog(f2t,P2t)
            hold on
            loglog(f2c,P2c)
            loglog(f2c,P2c*coeff_dTdC(i)^2)
            ylabel('PSD (S^2m^{-2}Hz^{-1})/(^{\circ}C^2Hz^{-1})')
            xlabel('f (Hz)')
            legend('dT','dC (pre calibration)','dC (post calibration)')
            title(sprintf('dT/dC=%.1f',coeff_dTdC(i)))
            grid on
            a3=subplot(3,1,3);
            plot(fctd.time,fctd.pressure,'k')
            hold on
            plot(fctd.time(i1),fctd.pressure(i1),'r','linewidth',2)
            ylabel('Pressure')
            datetick('x','keepticks','keeplimits')
            linkaxes([a1,a2],'x')
            pause
             clf
        end
    end
    
end
coeff_dTdC(coeff_dTdC>100)=NaN;
coeff_uC(coeff_uC>100)=NaN;

coeff_dTdC=nanmean(coeff_dTdC);
coeff_uC=nanmean(coeff_uC);
 else 
      coeff_dTdC=NaN;
      coeff_uC=NaN;
 end
end
