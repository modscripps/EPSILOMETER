 function mod_epsilometer_adjusted_chi(Meta_Data,Tnoise)

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


h_freq=Meta_Data.PROCESS.h_freq;
FPO7noise=Tnoise;

fc1=Meta_Data.PROCESS.fc1;
fc2=Meta_Data.PROCESS.fc2;

nfft=Meta_Data.PROCESS.nfft;
Fs_epsi=Meta_Data.PROCESS.Fs_epsi;

tscan=Meta_Data.PROCESS.tscan;
N_epsi=tscan.*Fs_epsi-mod(tscan*Fs_epsi,2);

dTdV(1)=Meta_Data.epsi.t1.dTdV; % define in mod_epsi_temperature_spectra
dTdV(2)=Meta_Data.epsi.t2.dTdV; % define in mod_epsi_temperature_spectra 


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
        Profile.chi2=zeros(nbscan,2).*nan;
        Profile.tg_fc2=zeros(nbscan,2).*nan;
        Profile.tg_flag2=zeros(nbscan,2).*nan;
        
        for s=1:nbscan % p is the scan index.
            if Profile.process_flag(s)==1
                scan.w=Profile.w(s);
                scan.ktemp=kt(Profile.s(s),Profile.t(s),Profile.pr(s));
                [~,indP] = sort(abs(Profile.P-Pr(s)));
                indP=indP(1);
                ind_Pr_epsi = find(Profile.epsitime<Profile.ctdtime(indP),1,'last');
                ind_scan = ind_Pr_epsi-N_epsi/2:ind_Pr_epsi+N_epsi/2; % ind_scan is even
                for c=1:nb_channels
                    wh_channel=channels{c};
                    switch wh_channel
                        case 't1'
                            scan.(wh_channel)=Profile.(wh_channel)(ind_scan).*dTdV(1); % time series in Celsius
                        case 't2'
                            scan.(wh_channel)=Profile.(wh_channel)(ind_scan).*dTdV(2); % time series in Celsius
                    end
                end
                
                [~,~,~,scan.chi(1),scan.tg_fc(1),~]=mod_efe_scan_chi(scan,'t1',Meta_Data,h_freq,FPO7noise);
                [~,~,~,scan.chi(2),scan.tg_fc(2),~]=mod_efe_scan_chi(scan,'t2',Meta_Data,h_freq,FPO7noise);
                

                Profile.chi(s,1)=scan.chi(1);
                Profile.chi(s,2)=scan.chi(2);
                Profile.tg_fc(s,1)=scan.tg_fc(1);
                Profile.tg_fc(s,2)=scan.tg_fc(2);
                Profile.tg_flag2(s,1)=1;
                Profile.tg_flag2(s,2)=1;


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





