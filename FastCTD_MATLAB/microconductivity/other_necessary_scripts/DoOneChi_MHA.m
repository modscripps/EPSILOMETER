function out=DoOneChi_MHA(FCTD,tlim,chi_param)
%function out=DoOneChi_MHA(FCTD,tlim,chi_param)
%Compute chi for the specified time period in tlim using the parameters in chi_param.  If the latter is not specified, defaults will be used.
%The FCTD structure must have had micro data processed by
%add_chi_microMHA_v2.m.

% BETHAN EDITS: 29/4/20
% - calculation of the power spectrum
% using data2 and datat seems to be causing problems. Only being used for
% plotting in this script and don't think it is therefore useful for the
% rest of the processing at the moment => edited the script so these
% calculations are only done when the length is greater than the length of
% the window and therefore the plotting is only done here too
% - for cases where data is too short, set out to NaN

if nargin < 3
    chi_param=FCTD_DefaultChiParam;
end

%Extract params from params struct
fs=chi_param.fs;
fslow=chi_param.fslow;
nsec=chi_param.nsec;
D=chi_param.D;
kmax=chi_param.kmax;
kmin=chi_param.kmin;
min_spd=chi_param.min_spd;
dt_sec=chi_param.dt_sec;
dTdC=chi_param.dTdC;
gain=chi_param.gain;
offset=chi_param.offset;
plotit=chi_param.plotit;

%Get indices
i1=find(FCTD.microtime>tlim(1)& FCTD.microtime<tlim(2));
i2=find(FCTD.time>tlim(1)&FCTD.time<tlim(2));

%Get indexed data
data=FCTD.ucon_corr(i1);
data2=FCTD.conductivity(i2); %This is low-pass C in
datat=FCTD.temperature(i2); %This is low-pass C in


%Mean fall rate
w=nanmean(FCTD.dPdt(i2));
spd=abs(w); %speed

%mean P
p=nanmean(FCTD.pressure(i2));

%High-freq
WINDOW=fs*nsec;

if length(data)<WINDOW || sum(~isnan(data))<2
    out.chi_stupid=NaN;
    out.chi1=NaN;
    out.w=NaN;
    out.pres=NaN;

else

    [P,f] = pwelch(diff(NANinterp(data))*fs,WINDOW,[],WINDOW,fs,'psd'); %Set NFFT equal to the window length

    %Compute chi from the standard formula as 6*D*(RMS dTdz)^2
    chi_stupid=6*D*dTdC^2*nanstd(diff(data)*fs/w).^2;
    chi1=nan;

    if spd > min_spd
        ichi=find(f/spd > kmin & f/spd < kmax);
        df=f(2)-f(1); %freq res
        %    chi1=6*D*nansum(P*spd*dTdC^2)*df/spd; %integrate wavenumber spectrum to get variance
        %11/2019 MHA: Caught two bugs: a missed factor of w, and also not integrating over a range.
        %P is the spectrum of dC/dt in (cond/s)^2/Hz.
        %chi1=6*D*dTdC^2*<dC/dz^2>
        %      =6*D*dTdC^2*<w^(-2)*dC/dt^2>
        %      =6*D*dTdC^2*w^(-2)*int(P)*df
        chi1=6*D*dTdC^2./spd^2*nansum(P(ichi))*df;
        %    chi1=6*D*nansum(P*spd*dTdC^2)*df/spd; %integrate wavenumber spectrum to get variance
    end

    %low-freq
    WINDOW=fslow*nsec;

    if length(data2)>WINDOW && length(datat)>WINDOW
        [P2,f2] = pwelch(diff(NANinterp(data2))*fslow,WINDOW,[],WINDOW,fslow,'psd'); %Set NFFT equal to window length

        %temperature
        [Pt,ft] = pwelch(diff(NANinterp(datat))*fslow,WINDOW,[],WINDOW,fslow,'psd'); %Set NFFT equal to window length

        %Plot
        if plotit
            %first difference correction, for plotting
            diffcorr=sinc(f2/fslow); %Compute first difference correction


            %Compute a batchelor spectrum or two for comparison
            eps1=1e-10;
            eps2=1e-8;
            eps3=1e-6;
            %     chi1=1e-10;
            %     chi2=1e-8;
            nu=1e-6;
            q=3.7;
            %     [k1,Pb1]=batchelor(eps1,chi_stupid,nu,D,q);
            %     [k2,Pb2]=batchelor(eps2,chi_stupid,nu,D,q);
            %     [ks,Pbs]=batchelor(eps3,chi_stupid,nu,D,q);
            [k1,Pb1]=batchelor(eps1,chi1,nu,D,q);
            [k2,Pb2]=batchelor(eps2,chi1,nu,D,q);
            [ks,Pbs]=batchelor(eps3,chi1,nu,D,q);

            figure(1)
            plot(FCTD.microtime(i1),data,FCTD.time(i2),data2)
            datetick
            freqline(tlim(1));
            freqline(tlim(2));

            figure(2)
            clf
            ax=MySubplot(.1,.1,0,.1,.1,0.1,1,3);
            axes(ax(1))
            plot(FCTD.time,FCTD.pressure,FCTD.time(i2),FCTD.pressure(i2),'r-','linewidth',2)
            axis ij
            %    xtloff
            datetick
            axes(ax(2))
            plot(FCTD.microtime,FCTD.ucon_corr,FCTD.microtime(i1),FCTD.ucon_corr(i1),'r-','linewidth',2)
            datetick
            axes(ax(3))
            plot(FCTD.time,FCTD.temperature,FCTD.time(i2),FCTD.temperature(i2),'r-','linewidth',2)
            datetick
            %    axes(ax(4))
            %    plot(FCTD.time,FCTD.salinity,FCTD.time(i2),FCTD.salinity(i2),'r-')
            %    datetick
            linkaxes(ax,'x')

            figure(3)
            clf
            loglog(f,P,f2,P2,f2,P2./diffcorr)

            ylabel('\Phi_{dCdz} / (S/m)^2/Hz')
            xlabel('\omega / Hz')
            %    ylim([1e-8 1e-4])
            %    xlim([1e-1 100])
            grid
            title(['w=' num2str(w) ' m/s'])
            shg

            %Now plot wavenumber spectrum of inferred temperature gradient.

            figure(4)
            clf
            loglog(f/spd,P*spd*dTdC^2,f2/spd,P2*spd*dTdC^2,ft/spd,Pt*spd,k1,Pb1,'k--',k2,Pb2,'k--',ks,Pbs,'k--')
            ylabel('\Phi_{dCdz} / (degC)/cpm')
            xlabel('m / cpm')
            legend('from micro','from C','from T','Batchelor')
            %    ylim([1e-7 1e-3])
            %    xlim([1e-1 100])
            grid
            %title(['TG wavenumber spec, w=' num2str(w) ' m/s, \chi_{stupid}=' num2str(chi_stupid)])
            title(['TG wavenumber spec, w=' num2str(w) ' m/s, \chi=' num2str(chi_stupid)])
            shg
            pause
        end
    end

    out.chi_stupid=chi_stupid;
    out.chi1=chi1;
    out.w=w;
    out.pres=p;

end
