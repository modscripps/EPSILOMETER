function [H1,H2]=create_velocity_tranfer_function(MS,min_epsi,Meta_Data)
% Create a transfer function using the Cross spectrum of u (from shear probes)
% and a3 or a1 ( it depends which axes influence most the shear probes).
% Turbulence_scans are raw spectra and time series from a single deployment
% with the lowest measured epsilon. It will serve to define the minimal noise floor 
% and the shape of the noise we can remove from every others scans. 
% The number of scans is undefined (user's discretion).
% 
%
% code inspired by Bethan Wynne-Cattanach's code
% written by Arnaud Le Boyer 02/06/2020.

listfile=dir(fullfile(Meta_Data.paths.profiles,'Turbulence_Profiles*.mat'));
listfilename=natsort({listfile.name});


Fs=Meta_Data.PROCESS.Fs_epsi;

average_min_eps1=min_epsi(1);
average_min_eps2=min_epsi(2);

id_min_scan1=cellfun(@(x) find(log10(x.epsilon(:,1))<average_min_eps1),MS,'un',0);
id_min_scan2=cellfun(@(x) find(log10(x.epsilon(:,2))<average_min_eps2),MS,'un',0);

count1=1;
count2=1;
scans1=[];
scans2=[];
for f=1:length(id_min_scan1)
    if (~isempty(id_min_scan1{f}))
        id_file=floor((f-1)./10)+1;
        load(fullfile(listfile(id_file).folder,listfilename{id_file}),sprintf('Profile%03i',f))
        eval(sprintf('Profile=Profile%03i;',f));
        for i=1:length(id_min_scan1{f})
            scans1{count1}.Pr=Profile.pr(id_min_scan1{f}(i));
            scans1{count1}.w=Profile.w(id_min_scan1{f}(i));
            scans1{count1} =mod_epsilometer_make_scan_v2(Profile,scans1{count1},Meta_Data);
            count1=count1+1;
        end
    end
    if (~isempty(id_min_scan2{f}))
        id_file=floor((f-1)./10)+1;
        load(fullfile(listfile(id_file).folder,listfilename{id_file}),sprintf('Profile%03i',f))
        eval(sprintf('Profile=Profile%03i;',f));
        for i=1:length(id_min_scan2{f})
            scans2{count2}.Pr=Profile.pr(id_min_scan2{f}(i));
            scans2{count2}.w=Profile.w(id_min_scan2{f}(i));
            scans2{count2} =mod_epsilometer_make_scan_v2(Profile,scans2{count2},Meta_Data);
            count2=count2+1;
        end
    end
end

%%
Nscans1=length(scans1);
Nscans2=length(scans2);
TF1=[];
TF2=[];
for n=1:Nscans1
    nfft=ceil(length(scans1{n}.s1)/3); %window length
    if mod(nfft,2)==0
        TF=zeros(Nscans1,ceil(nfft/2)+1);
    else
        TF=zeros(Nscans1,ceil(nfft/2));
    end
    
    u1 = scans1{n}.s1; %attention s1 is actually in m/s. Change was done in mod_som_make_scan_v2
    a3 = scans1{n}.a3;%attention a3 is actually in m/s^2. Change was done in mod_som_make_scan_v2
    [Pa3,fH] = pwelch(detrend(a3),nfft,[],nfft,Fs,'psd');
    [PCu1a3,~]=cpsd(detrend(u1),detrend(a3),nfft,[],nfft,Fs);
    TF1(n,:)=abs(PCu1a3)./Pa3;
end
for n=1:Nscans2
    nfft=ceil(length(scans2{n}.s2)/3); %window length
    if mod(nfft,2)==0
        TF=zeros(Nscans2,ceil(nfft/2)+1);
    else
        TF=zeros(Nscans2,ceil(nfft/2));
    end
    
    u2 = scans2{n}.s2; %attention s2 is actually in m/s. Change was done in mod_som_make_scan_v2
    a3 = scans2{n}.a3;%attention a3 is actually in m/s^2. Change was done in mod_som_make_scan_v2
    [Pa3,fH] = pwelch(detrend(a3),nfft,[],nfft,Fs,'psd');
    [PCu2a3,~]=cpsd(detrend(u2),detrend(a3),nfft,[],nfft,Fs);
    TF2(n,:)=abs(PCu2a3)./Pa3;

end

if ~isempty(TF1)
    H1=smoothdata(nanmean(TF1,1),'movmean',3);
    H1=interp1(fH,H1,Meta_Data.PROCESS.fe);
else
    H1=fH*0;
end
if ~isempty(TF2)
    H2=smoothdata(nanmean(TF2,1),'movmean',3);
    H2=interp1(fH,H2,Meta_Data.PROCESS.fe);
else
    H2=fH*0;
end



