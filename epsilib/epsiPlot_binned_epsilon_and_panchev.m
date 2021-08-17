function epsiPlot_binned_epsilon_and_panchev(Meta_Data,Profile,epsi_bin)

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
if nargin<3
    epsi_bin=10.^(-10:.5:-4.5);
end


% %% 
% listfile=dir(fullfile(Meta_Data.L1path,'Turbulence_Profiles*.mat'));
% listfilename=natsort({listfile.name});

% get channels
Fs_epsi=Meta_Data.AFE.FS;

tscan=Meta_Data.PROCESS.tscan;
N_epsi=tscan.*Fs_epsi-mod(tscan*Fs_epsi,2);

% Gravity  ... of the situation :)
G       = 9.81;
twoG= 2*G;
% calibration
Sv1=Meta_Data.AFE.s1.cal;
Sv2=Meta_Data.AFE.s2.cal;


%% create common k_axis
freq_min=0.1;    % 10 second^{-1}
freq_max=160;     % Nyquist=320 hz/2
average_speed=.5; %50 cm/s
k=(freq_min:freq_min:freq_max)./average_speed;


%% 
mkvis1=epsi_bin*0;
mspeed1=epsi_bin*0;
nbin1=epsi_bin*0;
mPsk1=zeros(length(epsi_bin),length(k));

mkvis2=epsi_bin*0;
mspeed2=epsi_bin*0;
nbin2=epsi_bin*0;
mPsk2=zeros(length(epsi_bin),length(k));


for i=1:length(epsi_bin)
        index1=find(Profile.epsilon(:,1)>=.5*epsi_bin(i) & Profile.epsilon(:,1)<1.5*epsi_bin(i));
        index2=find(Profile.epsilon(:,2)>=.5*epsi_bin(i) & Profile.epsilon(:,2)<1.5*epsi_bin(i));
%     index1=cellfun(@(x) find(x.epsilon(:,1)>=.5*epsi_bin(i) & x.epsilon(:,1)<1.5*epsi_bin(i)),Profile,'un',0);
%     index2=cellfun(@(x) find(x.epsilon(:,2)>=.5*epsi_bin(i) & x.epsilon(:,2)<1.5*epsi_bin(i)),Profile,'un',0);
%     notempty_index1=find(~isempty(index1));
%     notempty_index2=find(~isempty(index2));
    count_s1=1;
    count_s2=1;
    Psk1=0*k;
    kvis1=0;
    speed1=0;
    Psk2=0*k;
    kvis2=0;
    speed2=0;
    if ~isempty(notempty_index1)
        for id_profile=notempty_index1
           fprintf('epsi1 %1.2f, Pr=%i\n',log10(epsi_bin(i)),id_profile)
           if mod(id_profile,10)==0
                id_file=floor(id_profile/10);
            else
                id_file=floor(id_profile/10)+1;
            end
            load(fullfile(listfile(id_file).folder,listfilename{id_file}),sprintf('Profile%03i',id_profile))
            eval(sprintf('Profile=Profile%03i;',id_profile));
            scan.Cu1a.a3=Profile.Cu1a3;
            for id_scan=index1{id_profile}'
                Pr=Profile.pr(id_scan);
                scan.w=Profile.w(id_scan);
                scan.kvis=nu(Profile.s(id_scan),Profile.t(id_scan),Profile.pr(id_scan));
                [~,indP] = sort(abs(Profile.P-Pr));
                indP=indP(1);
                ind_Pr_epsi = find(Profile.epsitime<Profile.ctdtime(indP),1,'last');
                ind_scan = ind_Pr_epsi-N_epsi/2:ind_Pr_epsi+N_epsi/2; % ind_scan is even
                
%                 scan.a3=Profile.a3(ind_scan)*G; % time series in m.s^{-2}
                scan.s1=Profile.s1(ind_scan).*twoG./(Sv1.*scan.w); % time series in m.s^{-1}
                scan.s2=Profile.s1(ind_scan).*twoG./(Sv1.*scan.w); % time series in m.s^{-1}
                [~,~,~,P1,~,~,~,fe]=mod_efe_scan_epsilon(scan,'s1','a3',Meta_Data);
                
                Psk1=Psk1+interp1(fe/scan.w,P1,k);
                kvis1=kvis1+scan.kvis;
                speed1=speed1+scan.w;
                count_s1=count_s1+1;
            end
        end
        nbin1(i)=count_s1;
        mPsk1(i,:)=Psk1./count_s1;
        mkvis1(i)=kvis1./count_s1;
        mspeed1(i)=speed1./count_s1;
    end
    if ~isempty(notempty_index2)
        for id_profile=notempty_index2
            fprintf('epsi2 %1.2f, Pr=%i\n',log10(epsi_bin(i)),id_profile)
            if mod(id_profile,10)==0
                id_file=floor(id_profile/10);
            else
                id_file=floor(id_profile/10)+1;
            end
            load(fullfile(listfile(id_file).folder,listfilename{id_file}),sprintf('Profile%03i',id_profile))
            eval(sprintf('Profile=Profile%03i;',id_profile));
            scan.Cu2a.a3=Profile.Cu2a3;
            for id_scan=index2{id_profile}'
                Pr=Profile.pr(id_scan);
                scan.w=Profile.w(id_scan);
                scan.kvis=nu(Profile.s(id_scan),Profile.t(id_scan),Profile.pr(id_scan));
                [~,indP] = sort(abs(Profile.P-Pr));
                indP=indP(1);
                ind_Pr_epsi = find(Profile.epsitime<Profile.ctdtime(indP),1,'last');
                ind_scan = ind_Pr_epsi-N_epsi/2:ind_Pr_epsi+N_epsi/2; % ind_scan is even
                
%                 scan.a3=Profile.a3(ind_scan)*G; % time series in m.s^{-2}
                % this is bad but I am too tired to fixed now
                % I copy s1 and s2 so mod_efe_scan_acceleration can run.
                % Because I just changed it and I did not think about
                % binned_epsi.
                scan.s1=Profile.s2(ind_scan).*twoG./(Sv2.*scan.w); % time series in m.s^{-1}
                scan.s2=Profile.s2(ind_scan).*twoG./(Sv2.*scan.w); % time series in m.s^{-1}
%                 wh_channel='a3';
%                 [~,scan.Cu2a.(wh_channel),~,~,~]=mod_efe_scan_coherence(scan,wh_channel,Meta_Data);
                
                [~,~,~,P2,~,~,~,fe]=mod_efe_scan_epsilon(scan,'s2','a3',Meta_Data);
                
                Psk2=Psk2+interp1(fe/scan.w,P2,k);
                kvis2=kvis2+scan.kvis;
                speed2=speed2+scan.w;
                count_s2=count_s2+1;
            end
        end
        nbin2(i)=count_s2;
        mPsk2(i,:)=Psk2./count_s2;
        mkvis2(i)=kvis2./count_s2;
        mspeed2(i)=speed2./count_s2;
    end
end


%% plotting session

mspeed1(mspeed1==0)=nan;
mspeed2(mspeed2==0)=nan;
mkvis1(mkvis1==0)=nan;
mkvis2(mkvis2==0)=nan;
w_th1=nanmean(mspeed1);
w_th2=nanmean(mspeed2);
kvis_th1=nanmean(mkvis1);
kvis_th2=nanmean(mkvis2);

% get shear channel average noise to compute epsi
logf=log10(fe);
h_freq=get_filters_MADRE(Meta_Data,fe);
shearnoise=load(fullfile(Meta_Data.CALIpath,'shear_noise.mat'),'n0s','n1s','n2s','n3s');
n0s=shearnoise.n0s; n1s=shearnoise.n1s; n2s=shearnoise.n2s; n3s=shearnoise.n3s;
snoise=10.^(n0s+n1s.*logf+n2s.*logf.^2+n3s.*logf.^3);
TFshear1=(Sv1.*w_th1/twoG).^2 .* h_freq.shear.* haf_oakey(fe,w_th1);
TFshear2=(Sv2.*w_th1/twoG).^2 .* h_freq.shear.* haf_oakey(fe,w_th2);
snoise_k1= (2*pi*fe./w_th1).^2 .* snoise.*w_th1./TFshear1;        % T1_k spec  as function of k
snoise_k2= (2*pi*fe./w_th2).^2 .* snoise.*w_th2./TFshear2;        % T1_k spec  as function of k


mkvis1(isnan(mkvis1))=nanmean(mkvis1);
mkvis2(isnan(mkvis2))=nanmean(mkvis2);

for i=1:length(epsi_bin)
    [kpan,Ppan]=panchev(epsi_bin(i),mkvis1(i));
    Ppan1(i,:)=interp1(kpan,Ppan,k);
    [kpan,Ppan]=panchev(epsi_bin(i),mkvis2(i));
    Ppan2(i,:)=interp1(kpan,Ppan,k);
end



F1=figure;
set(F1,'position',[100 50 1300 700])
ax1=axes('position',[.1 .1 .55 .82]);
cmap = parula(length(epsi_bin)+1);

pp = loglog(ax1,k,Ppan1,'color',[1 1 1]*.7,'linewidth',1);
ll=pp;

mink=k(find(~isnan(nansum(mPsk1,1)),1,'first'));
maxk=k(find(~isnan(nansum(mPsk1,1)),1,'last'));

mPsk1(mPsk1==0) = nan
minP=min(mPsk1(:));
maxP=max(mPsk1(:));

hold on
for i=1:length(epsi_bin)
    indk=find(~isnan(Ppan1(i,:)),1,'first');
    if k(indk)<=mink
        text(mink, ...
            Ppan1(i,indk),...
            sprintf('10^{%1.1f}',log10(epsi_bin(i))),...
            'fontsize',20,'Parent',ax1)
    else
        text(k(indk), ...
            Ppan1(i,indk),...
            sprintf('10^{%1.1f}',log10(epsi_bin(i))),...
            'fontsize',20,'Parent',ax1)
    end
    ll(i)=loglog(ax1,k,mPsk1(i,:),'Color',cmap(i,:),'linewidth',2);
    loglog(ax1,fe./w_th1,10*snoise_k1,'Color',[.2 .2 .2],'linewidth',2);
    
end

xlim(ax1,[mink maxk])
ylim(ax1,[minP maxP])
grid(ax1,'on')
xlabel(ax1,'k [cpm]')
ylabel(ax1,'[s^{-2} / cpm]')
set(ax1,'fontsize',20)
title(ax1,sprintf('%s - %s - shear1',Meta_Data.mission,Meta_Data.deployment))

nanind=find(nansum(mPsk1,2)>0);
legend_string=arrayfun(@(x) sprintf('log_{10}(\\epsilon)~%2.1f',x),log10(epsi_bin),'un',0);
hl=gridLegend(ll(nanind),2,legend_string{nanind},'location','northwest');
set(hl,'fontsize',5)
hl.Position=[.01 .75 .25 .25];

ax2=axes('position',[.7 .1 .25 .82]);

plot(ax2,log10(epsi_bin),nbin1,'-+k')
set(ax2,'fontsize',20)
xlabel(ax2,'log_{10} \epsilon')
xlim([min(log10(epsi_bin)) max(log10(epsi_bin))])
title(ax2,'PDF \epsilon_1')
print(F1,'-dpng2',fullfile(Meta_Data.L1path,sprintf('%s_%s_shear1.png',Meta_Data.mission,Meta_Data.deployment)))

%% EPSILON 2
F2=figure;
set(F2,'position',[100 50 1300 700])
ax1=axes('position',[.1 .1 .55 .82]);
cmap = parula(length(epsi_bin)+1);

pp = loglog(ax1,k,Ppan2,'color',[1 1 1]*.7,'linewidth',1);
ll=pp;

mink=k(find(~isnan(nansum(mPsk2,1)),1,'first'));
maxk=k(find(~isnan(nansum(mPsk2,1)),1,'last'));

mPsk2(mPsk2==0) = nan
minP =min(mPsk2(:));
maxP=max(mPsk2(:));

hold on
for i=1:length(epsi_bin)
    indk=find(~isnan(Ppan2(i,:)),1,'first');
    if k(indk)<=mink
        text(mink, ...
            Ppan2(i,indk),...
            sprintf('10^{%1.1f}',log10(epsi_bin(i))),...
            'fontsize',20,'Parent',ax1)
    else
        text(k(indk), ...
            Ppan2(i,indk),...
            sprintf('10^{%1.1f}',log10(epsi_bin(i))),...
            'fontsize',20,'Parent',ax1)
    end
    ll(i)=loglog(ax1,k,mPsk2(i,:),'Color',cmap(i,:),'linewidth',2);
    loglog(ax1,fe./w_th2,10*snoise_k2,'Color',[.2 .2 .2],'linewidth',2);
    
end

xlim(ax1,[mink maxk])
ylim(ax1,[minP maxP])
grid(ax1,'on')
xlabel(ax1,'k [cpm]')
ylabel(ax1,'[s^{-2} / cpm]')
set(ax1,'fontsize',20)
title(ax1,sprintf('%s - %s - shear2',Meta_Data.mission,Meta_Data.deployment))

nanind=find(nansum(mPsk1,2)>0);
legend_string=arrayfun(@(x) sprintf('log_{10}(\\epsilon)~%2.1f',x),log10(epsi_bin),'un',0);
hl=gridLegend(ll(nanind),2,legend_string{nanind},'location','northwest');
set(hl,'fontsize',5)
hl.Position=[.01 .75 .25 .25];

ax2=axes('position',[.7 .1 .25 .82]);

plot(ax2,log10(epsi_bin),nbin2,'-+k')
set(ax2,'fontsize',20)
xlabel(ax2,'log_{10} \epsilon')
xlim([min(log10(epsi_bin)) max(log10(epsi_bin))])
title(ax2,'PDF \epsilon_2')
print(F2,'-dpng2',fullfile(Meta_Data.L1path,sprintf('%s_%s_shear2.png',Meta_Data.mission,Meta_Data.deployment)))


