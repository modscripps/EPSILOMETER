function mod_epsilometer_adjusted_chi(Meta_Data,Tnoise1,Tnoise2)

%  Profile structure for Micro Structure. Inside Profile you ll find
%  temperature spectra in degC Hz^-1
%  Horizontal  velocity spectra in m^2/s^-2 Hz^-1
%  Acceleration/speed spectra in s^-1 Hz^-1 
%
%  Created by Arnaud Le Boyer on 7/28/18.
listfile=dir(fullfile(Meta_Data.paths.profiles,'Turbulence_Profiles*.mat'));
listfilename=natsort({listfile.name});

%% get channels
channels=Meta_Data.PROCESS.channels;
nb_channels=length(channels);

inda3=find(cellfun(@(x) strcmp(x,'a3'),channels));
inds1=find(cellfun(@(x) strcmp(x,'s1'),channels));
inds2=find(cellfun(@(x) strcmp(x,'s2'),channels));


fc1=Meta_Data.PROCESS.fc1;
fc2=Meta_Data.PROCESS.fc2;

nfft=Meta_Data.PROCESS.nfft;
Fs_epsi=Meta_Data.PROCESS.Fs_epsi;

tscan=Meta_Data.PROCESS.tscan;
N_epsi=tscan.*Fs_epsi-mod(tscan*Fs_epsi,2);

%% Gravity  ... of the situation :)
G       = 9.81;
twoG= 2*G;
Sv1=Meta_Data.epsi.s1.Sv;
Sv2=Meta_Data.epsi.s2.Sv;



%  fc1 and fc2 define the frequency range integration for qc.
count=1;
sav_var_name=[];
% TODO do not forget to change the listfile for loop back to length(listfile)
for f=1:length(listfile)
    load(fullfile(listfile(f).folder,listfilename{f}),'nb_profile_perfile')
    for p=1:nb_profile_perfile
        load(fullfile(listfile(f).folder,listfilename{f}),sprintf('Profile%03i',count))
        fprintf('%s:%s\r\n',fullfile(listfile(f).folder,listfilename{f}),sprintf('Profile%03i',count));
        eval(sprintf('Profile=Profile%03i;',count));
        
        Pr=Profile.pr;
        nbscan=Profile.nbscan;
        
        % define a Pressure axis to an which I will compute epsilon and chi.
        %  The spectra will be nfft long centered around P(z) +/- tscan/2.
        %
        %initialize process flags
        Profile.sh_qcflag=zeros(nbscan,2).*nan;
        Profile.epsilonTF=zeros(nbscan,2).*nan;
        Profile.sh_fcTF=zeros(nbscan,2).*nan;

%         TFnoise=@(x,y) (interp1(fH,y,x));
        scan.Cu1a.a1=Profile.Cu1a1;
        scan.Cu1a.a2=Profile.Cu1a2;
        scan.Cu1a.a3=Profile.Cu1a3;
        scan.Cu2a.a1=Profile.Cu2a1;
        scan.Cu2a.a2=Profile.Cu2a2;
        scan.Cu2a.a3=Profile.Cu2a3;

        for s=1:nbscan % p is the scan index.
            if Profile.process_flag(s)==1
                scan.w=Profile.w(s);
                scan.kvis=nu(Profile.s(s),Profile.t(s),Profile.pr(s));
                [~,indP] = sort(abs(Profile.P-Pr(s)));
                indP=indP(1);
                ind_Pr_epsi = find(Profile.epsitime<Profile.ctdtime(indP),1,'last');
                ind_scan = ind_Pr_epsi-N_epsi/2:ind_Pr_epsi+N_epsi/2; % ind_scan is even
                
                wh_channels=channels{inda3};
                
                scan.(wh_channels)=Profile.(wh_channels)(ind_scan)*G; % time series in m.s^{-2}
                [scan.P.(wh_channels),~] = pwelch(detrend(scan.(wh_channels)),nfft,[],nfft,Fs_epsi,'psd');
                
                u1_vibration=scan.P.(wh_channels).*H1;
                u2_vibration=scan.P.(wh_channels).*H2;
                
                scan.s1=Profile.s1(ind_scan).*twoG./(Sv1.*scan.w); % time series in m.s^{-1}
                [P1,~,~,~,scan.epsilon(1),scan.sh_fc(1),~]=mod_efe_scan_epsilon_withTF(scan,u1_vibration,'s1',Meta_Data);
                
                scan.s2=Profile.s2(ind_scan).*twoG ./(Sv2.*scan.w); % time series in m.s^{-1}
                [P2,~,~,~,scan.epsilon(2),scan.sh_fc(2),fe]=mod_efe_scan_epsilon_withTF(scan,u2_vibration,'s2',Meta_Data);
                
                % get the ratio of U_obs / U_vibration (from low epsilon region)
                qc_flag1= log10(P1./ smoothdata(u1_vibration,'movmean',5));
                qc_flag2= log10(P2./ smoothdata(u2_vibration,'movmean',5));
                
                qc_flag1=nanmean(qc_flag1(fe>fc1 & fe<fc2));
                qc_flag2=nanmean(qc_flag2(fe>fc1 & fe<fc2));
                
                Profile.sh_qcflag(s,:)=[qc_flag1 qc_flag2];

                Profile.epsilonTF(s,1)=scan.epsilon(1);
                Profile.epsilonTF(s,2)=scan.epsilon(2);
                Profile.sh_fcTF(s,1)=scan.sh_fc(1);
                Profile.sh_fcTF(s,2)=scan.sh_fc(2);


            end
        end
        eval(sprintf('Profile%03i=Profile;',count));
        clear Profile;
        sav_var_name=[sav_var_name sprintf(',''Profile%03i''',count)];
        count=count+1;
    end
    save_file=fullfile(Meta_Data.paths.profiles, ...
        ['Turbulence_Profiles' num2str((f-1)*10) '.mat']);
    cmd=['save(''' save_file '''' sav_var_name ',''nb_profile_perfile'')'];
    eval(cmd);
    sav_var_name=[];
end





