function [Coh_s1a,Coh_s2a,sumCoh_s1,sumCoh_s2,fe]=mod_efe_scan_multivariate_coherence(scan,nb_accel_channel,nfft,dof_coh,Fs)
% Correcting the shear channels using a multivariate acceleration
% correction: I use all the acceleration channnels to get a combined
% coherence that removed from the shear channel. 
%
% Input: scan, meta_data  
%        
% Compute the coherence over the whole profile using all 3 axis (we could add more) 
% over the 1./tsan:Fs frequency frequency axis with nfft samples.
%
% I changed the design of the correction. 
% Now I compute the multivariate coherence along the whole profile.
% The coherence is computed using the dof_coh in meta_data. If we follow
% rolf lueck recommandation dof_coh=15 is good. 
% This means that we need to "slide" our multivariate approach along the
% profile using segments with a length = dof_coh * NFFT and obtain the
% multivariate coherence for a freqeuncy array = Meta_Data.PROCESS.fe

% written by aleboyer@ucsd.edu 10/01/2021

% when the probes are oriented down/upward
% a3 should be alingned with the shear
% a1 is the other horizontal accell
% a2 is the vertical 
list_accel={'a3_g','a1_g','a2_g'};

Csi=zeros(nb_accel_channel,nfft/2+1);
Csj=zeros(nb_accel_channel,nfft/2+1);
Aij=zeros(nb_accel_channel,nb_accel_channel,nfft/2+1);
MC=zeros(nb_accel_channel,nb_accel_channel,nfft/2+1);

df=Fs./nfft;
if isfinite(scan.s1_volt)
    % lets do this 
    [Pv1,fe]=pwelch(detrend(scan.s1_volt),nfft,[],nfft,Fs);

    for i=1:nb_accel_channel% nb of accel channel
        wh_accel_i=list_accel{i};
        accel_i=detrend(scan.(wh_accel_i));
        Csi(i,:)=cpsd(detrend(scan.s1_volt), accel_i,...
                      nfft,[],nfft,Fs);

        for j=1:nb_accel_channel %nb of accel channel
            wh_accel_j=list_accel{i};
            accel_j=detrend(scan.(wh_accel_j));

            Csj(j,:)=cpsd(detrend(scan.s1_volt), accel_j,...
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
    Coh=squeeze(nansum(nansum(MC,1),2))*df;      
    Pv1_co=Pv1-Coh;
  
end

if isfinite(scan.s2_volt)
    % lets do this 
    [Pv2,~]=pwelch(detrend(scan.s2_volt),nfft,[],nfft,Fs);

    for i=1:nb_accel_channel% nb of accel channel
        wh_accel_i=list_accel{i};
        accel_i=detrend(scan.(wh_accel_i));
        Csi(i,:)=cpsd(detrend(scan.s2_volt), accel_i,...
                      nfft,[],nfft,Fs);

        for j=1:nb_accel_channel %nb of accel channel
            wh_accel_j=list_accel{i};
            accel_j=detrend(scan.(wh_accel_j));

            Csj(j,:)=cpsd(detrend(scan.s2_volt), accel_j,...
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
    Coh=squeeze(nansum(nansum(MC,1),2))*df;      
    Pv2_co=Pv2-Coh;
end