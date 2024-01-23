function [Coh_s1a,Coh_s2a]=mod_efe_profile_multivariate_coherence(Profile,Pr,Meta_Data)

% Correcting the shear channels using a multivariate acceleration
% correction: I use all the acceleration channnels to get a combined
% coherence that removed from the shear channel. 
%
% Input: Profile, meta_data  
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

nfft=Meta_Data.PROCESS.nfft;
Fs=Meta_Data.PROCESS.Fs_epsi;
dof_coh=Meta_Data.PROCESS.dof_coh;
nb_accel_channel=numel(list_accel);

Coh_s1a=zeros(numel(Pr),nfft/2+1);
Coh_s2a=zeros(numel(Pr),nfft/2+1);
N_epsi      = (dof_coh-1)*Meta_Data.PROCESS.nfft;

for p=1:numel(Pr)
    fprintf("Multivariate Pr=%3.1f, Pr(end)=%3.1f \r\n",Pr(p),Pr(end))
    % Find indices of ctd and epsi data that make up this scan
    [~,indP]    = sort(abs(Profile.ctd.P-Pr(p)));
    indP        = indP(1);
    ind_Pr_epsi = find(Profile.time_s<Profile.ctd.time_s(indP),1,'last');
    ind_scan    = ind_Pr_epsi-N_epsi/2:ind_Pr_epsi+N_epsi/2-1; % ind_scan is even
    % ALB case where we are to close to the beginning of the profile 
    if ind_scan(1)<=0
        ind_scan=ind_scan-ind_scan(1)+1;
    end
    % ALB case where we are to close to the end of the profile
    if ind_scan(end)>=numel(Profile.time_s)
        ind_scan=ind_scan-(ind_scan(end)-numel(Profile.time_s));
    end
    
    Csi=zeros(nb_accel_channel,nfft/2+1);
    Csj=zeros(nb_accel_channel,nfft/2+1);
    Aij=zeros(nb_accel_channel,nb_accel_channel,nfft/2+1);
    MC=zeros(nb_accel_channel,nfft/2+1);
    
    df=Fs./nfft;
    if isfinite(Profile.s1_volt)
        for i=1:nb_accel_channel% nb of accel channel
            wh_accel_i=list_accel{i};
            det_accel_i=detrend(Profile.(wh_accel_i)(ind_scan));
            Csi(i,:)=cpsd(detrend(Profile.s1_volt(ind_scan)), det_accel_i, ...
                nfft,[],nfft,Fs);
            Csj(i,:)=cpsd(det_accel_i,detrend(Profile.s1_volt(ind_scan)),  ...
                    nfft,[],nfft,Fs);
            
            for j=1:nb_accel_channel %nb of accel channel
                wh_accel_j=list_accel{j};
                det_accel_j=detrend(Profile.(wh_accel_j)(ind_scan));
                
                
                Aij(i,j,:)=cpsd(det_accel_i,det_accel_j, ...
                    nfft,[],nfft,Fs);
            end
        end
        
        for i=1:nb_accel_channel% nb of accel channel
            for j=1:nb_accel_channel
                MC(i,:)=MC(i,:)+Csi(j,:)./ squeeze(Aij(j,i,:)).'.*conj(Csj(j,:));
            end
        end
        % sum the MC (Mutlivariate Coeficients to get the correction)
        % Not sure about the df. It is in the Goodman2006 paper but I have a
        % doubt wether it is already included in the output of cpsd.
        % With the df It seems to do a good job in the correction.
        Coh_s1a(p,:)=squeeze(nansum(MC))*df;
    end
    
    if isfinite(Profile.s2_volt)
        
        Csi=zeros(nb_accel_channel,nfft/2+1);
        Csj=zeros(nb_accel_channel,nfft/2+1);
        Aij=zeros(nb_accel_channel,nb_accel_channel,nfft/2+1);
        MC=zeros(nb_accel_channel,nfft/2+1);

        for i=1:nb_accel_channel% nb of accel channel
            wh_accel_i=list_accel{i};
            det_accel_i=detrend(Profile.(wh_accel_i)(ind_scan));
            Csi(i,:)=cpsd(detrend(Profile.s2_volt(ind_scan)), det_accel_i,...
                nfft,[],nfft,Fs);
            Csj(i,:)=cpsd(det_accel_i,detrend(Profile.s2_volt(ind_scan)),  ...
                nfft,[],nfft,Fs);
            
            for j=1:nb_accel_channel %nb of accel channel
                wh_accel_j=list_accel{i};
                det_accel_j=detrend(Profile.(wh_accel_j)(ind_scan));
                
                Aij(i,j,:)=cpsd(det_accel_i,det_accel_j, ...
                    nfft,[],nfft,Fs);
            end
        end
        for i=1:nb_accel_channel% nb of accel channel
            for j=1:nb_accel_channel
                MC(i,:)=MC(i,:)+Csi(j,:)./ squeeze(Aij(j,i,:)).'.*conj(Csj(j,:));
            end
        end
        
        % sum the MC (Mutlivariate Coeficients to get the correction)
        % Not sure about the df. It is in the Goodman2006 paper but I have a
        % doubt wether it is already included in the output of cpsd.
        % With the df It seems to do a good job in the correction.
        Coh_s2a(p,:)=squeeze(nansum(MC))*df;

    end
%     check_multivariate
end

%ALB following the SCOR group recommendation average 20 scan
dof_coherence=20;
Coh_s2a=abs(smoothdata(Coh_s2a,'movmedian',dof_coherence));
Coh_s1a=abs(smoothdata(Coh_s1a,'movmedian',dof_coherence));


disp("multivariate done.")
