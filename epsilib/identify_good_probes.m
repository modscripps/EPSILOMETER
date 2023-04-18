root='/Volumes/GoogleDrive/Shared drives/MOD/MOD Cruises & Experiments/TLC/Cruises/Pilot_2023/05_processed_data/';


list_dep=dir([root,'*epsi*']);
%%
nb_profile=0;
for d=1:length(list_dep)
    dep_name=fullfile(list_dep(d).folder,list_dep(d).name);
    listprofile=dir(fullfile(dep_name,'profiles','Profile*.mat'));
    for f=1:length(listprofile)
        nb_profile=nb_profile+1;
        filename=fullfile(listprofile(f).folder,listprofile(f).name);
        disp(filename)
        load(filename)
        dep_series(nb_profile)=d;
        prof_series(nb_profile)=f;
        SN_shear(nb_profile,1)=str2num(Profile.Meta_Data.AFE.s1.SN);
        SN_shear(nb_profile,2)=str2num(Profile.Meta_Data.AFE.s2.SN);
        SN_temp(nb_profile,1)=str2num(Profile.Meta_Data.AFE.t1.SN);
        SN_temp(nb_profile,2)=str2num(Profile.Meta_Data.AFE.t2.SN);

        meanepsi(nb_profile,1)=nanmean(log10(Profile.epsilon(:,1)));
        meanepsi(nb_profile,2)=nanmean(log10(Profile.epsilon(:,2)));
        
        meantemp(nb_profile,1)=nanmean(log10(Profile.chi(:,1)));
        meantemp(nb_profile,2)=nanmean(log10(Profile.chi(:,2)));

    end

end

%%
close all
figure(1)
ax(1)=subplot(411);
plot(SN_shear(:,1),'b','linewidth',2)
hold on
plot(SN_shear(:,2),'r','linewidth',2)
plot(SN_temp(:,1),'k','linewidth',2)
plot(SN_temp(:,2),'g','linewidth',2)
legend(ax(1),'s1','s2','t1','t2')
ax(2)=subplot(412);
plot(meanepsi(:,1),'b','linewidth',2)
hold on
plot(meanepsi(:,2),'r','linewidth',2)
plot(meantemp(:,1),'k','linewidth',2)
plot(meantemp(:,2),'g','linewidth',2)
legend(ax(2),'s1','s2','t1','t2')

ax(3)=subplot(413);
plot(dep_series,'linewidth',2)

ax(4)=subplot(414);
plot(prof_series,'linewidth',2)

linkaxes(ax,'x')

%%
d=2;
f=39;
dep_name=fullfile(list_dep(d).folder,list_dep(d).name);
listprofile=dir(fullfile(dep_name,'profiles','Profile*.mat'));

filename=fullfile(listprofile(f).folder,listprofile(f).name);
load(filename)
figure(2)
ax(1)=subplot(121);
semilogx(Profile.epsilon_co,Profile.pr);axis ij
title('shear')
ax(2)=subplot(122);
semilogx(Profile.chi,Profile.pr);axis ij
title('temp')
