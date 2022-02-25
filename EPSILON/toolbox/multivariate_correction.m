function [Pv,Pv_co,fe]=multivariate_correction(signal2correct,vibration,nfft,Fs)
% Correcting the shear channels using a multivariate acceleration
% correction: I use all the acceleration channnels to get a combined
% coherence that removed from the shear channel. 
%
% Input: signaltocorrect: 1D timeseries of the signal you want to correct
%        vibration      : 2D (timestamp,channels) timseries of "vibrations"
%                         you want to use to correct signaltocorrect  
%        
% Compute the coherence over the whole profile using all 3 axis (we could add more) 
% over the 1./tsan:Fs frequency frequency axis with nfft samples.

% written by aleboyer@ucsd.edu 10/01/2021

% when the probes are oriented down/upward
% a3 should be alingned with the shear
% a1 is the other horizontal accell
% a2 is the vertical 

if size(vibration,2)>size(vibration,1)
    vibration=vibration.';
end

% detrend signal to correct before fft
detrend_signal2correct=detrend(signal2correct);

% detrend vibration timeseries before fft
detrend_vibration=detrend(vibration);


nb_vibration_channel=size(vibration,2);
L_signal=length(signal2correct);
dof=2*L_signal/nfft-1; % for 50% overlap segments

Csi=zeros(nb_vibration_channel,nfft/2+1);
Csj=zeros(nb_vibration_channel,nfft/2+1);
Aij=zeros(nb_vibration_channel,nb_vibration_channel,nfft/2+1);
MC=zeros(nb_vibration_channel,nb_vibration_channel,nfft/2+1);

df=Fs./nfft;
    % lets do this 
    [Pv,fe]=pwelch(detrend_signal2correct,nfft,[],nfft,Fs);

    for i=1:nb_vibration_channel% nb of accel channel
        accel_i=detrend_vibration(:,i);
        Csi(i,:)=cpsd(detrend_signal2correct, accel_i,...
                      nfft,[],nfft,Fs);

        for j=1:nb_vibration_channel %nb of accel channel
            accel_j=detrend_vibration(:,j);

            Csj(j,:)=cpsd(detrend_signal2correct, accel_j,...
                nfft,[],nfft,Fs);

            Aij(i,j,:)=cpsd(accel_i,accel_j, ...
                            nfft,[],nfft,Fs);
                        
            MC(i,j,:)=Csi(i,:).*conj(Csj(j,:))./ ...
                      squeeze(Aij(i,j,:)).';
            
        end  
    end
    % sum the MC (Mutlivariate Coeficients to get the correction)
    % Not sure about the df. It is in the Goodman2006 paper but I have a
    % doubt wether it is already included in the output of cpsd. 
    % With the df It seems to do a good job in the correction.
    Coh=squeeze(nansum(nansum(MC,1),2));      
    bias=1-1.02*nb_vibration_channel./dof;
    Pv_co=(Pv-Coh)./bias;
  
end
