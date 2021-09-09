
function epsiPlot_binned_epsilon_and_panchev(Profile,epsi_bin)

%  compute averaged shear spectra classed by epsilon values
%
%  input: 
% . MS: structure with the chi values and the temperature gradiant spectra
% . epsi_bin: bin of epsilon values
%
%  Created by Arnaud Le Boyer on 7/28/18.
%  Aug. 2021 - Nicole Couto adapted mod_epsilometer_binned_epsilon to work
%  with epsi_class and Profile structures

%% create epsi_bin or include min and max epsi value to epsi_bin
if nargin<2
    epsi_lim=10.^(-10:.5:-4.5);
    epsi_bin = [epsi_lim(1:end-1);epsi_lim(2:end)].';
end

%% create common k_axis
freq_min=0.1;    % 10 second^{-1}
freq_max=160;     % Nyquist=320 hz/2
average_speed=.5; %50 cm/s
k=(freq_min:freq_min:freq_max)./average_speed;

for b=1:length(epsi_bin)
    
    Ps1_interp = [];
    Ps2_interp = [];
    
    % Find values of epsilon within this bin
    inBin1{b} = Profile.epsilon_co(:,1)>=epsi_bin(b,1) & Profile.epsilon_co(:,1)<=epsi_bin(b,2);
    inBin2{b} = Profile.epsilon_co(:,2)>=epsi_bin(b,1) & Profile.epsilon_co(:,1)<=epsi_bin(b,2);
    
    % Find average kinematic viscosity of these data
    kvis1 = nanmean(Profile.kvis(inBin1{b}));
    kvis2 = nanmean(Profile.kvis(inBin2{b}));
    
    % Find Panchev curve
    [~,panchev1{b}] = panchev(nanmean(epsi_bin(b,:)),kvis1,k);
    [~,panchev2{b}] = panchev(nanmean(epsi_bin(b,:)),kvis2,k);
    
    % Find average shear spectrum - need to interpolate Profile.k (which is
    % different for each scan) onto k
    Ps1 = Profile.Ps_shear_co_k.s1(inBin1{b},:);
    k1 = Profile.k(inBin1{b},:);
    Ps2 = Profile.Ps_shear_co_k.s2(inBin2{b},:);
    k2 = Profile.k(inBin2{b},:);
    for ii=1:sum(inBin1{b})
        Ps1_interp(ii,:) = interp1(k1(ii,:),Ps1(ii,:),k);
    end
    for ii=1:sum(inBin2{b})
        Ps2_interp(ii,:) = interp1(k2(ii,:),Ps2(ii,:),k);
    end
    
    if ~isempty(Ps1_interp)
        Psh_avg1{b} = nanmean(Ps1_interp);
    else
        Psh_avg1{b} = nan(size(k));
    end
    if ~isempty(Ps2_interp)
        Psh_avg2{b} = nanmean(Ps2_interp);
    else
        Psh_avg2{b} = nan(size(k));
    end
      
end

%% Make a figure for each probe
cols = jet(length(epsi_bin));

figure('position',[0,0,19.8,13.2])

subtightplot(2,1,1)
for b=1:length(epsi_bin)
    legendString = sprintf('%2.1E - %2.1E, %3.0f',epsi_bin(b,1), epsi_bin(b,2),sum(inBin1{b}));
    p(b) = loglog(k,panchev1{b},'linewidth',4,'color',cols(b,:),'displayname',legendString);
    hold on
end
for b=1:length(epsi_bin)
    loglog(k,Psh_avg1{b},'linewidth',2,'color',cols(b,:));
end
legend(p,'location','eastoutside')
title('shear channel 1')
set(gca,'ylim',[1e-6,1e-1])

subtightplot(2,1,2)
for b=1:length(epsi_bin)
    legendString = sprintf('%2.1E - %2.1E, %3.0f',epsi_bin(b,1), epsi_bin(b,2),sum(inBin2{b}));
    p(b) = loglog(k,panchev2{b},'linewidth',4,'color',cols(b,:),'displayname',legendString);
    hold on
end
for b=1:length(epsi_bin)
    loglog(k,Psh_avg2{b},'linewidth',2,'color',cols(b,:));
end
legend(p,'location','eastoutside')
title('shear channel 2')
set(gca,'ylim',[1e-6,1e-1])



