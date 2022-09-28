% function concat_MS(path)
%concatanate Turbulence_Profiles.mat from the path 

path='/Volumes/GoogleDrive/Shared drives/MOD data 2019 MISO-BOB 2/epsi_data/mbob/epsifish_1/d1/L1/';
listfile=dir(fullfile(path,'Turbulence_Profile*mat'));
MS1=[];
for f=1:length(listfile)
    filepath=fullfile(listfile(f).folder,listfile(f).name);
    load(filepath)
    MS1=[MS1 MS];
end

%%
load /Volumes/GoogleDrive/'Shared drives'/'MOD data 2019 MISO-BOB 2'/epsi_data/mbob/epsifish_1/d1/L1/Profiles_d1.mat
%%
epsi1=[];
epsi2=[];
Pr_id=[];
scan_id=[];
for f=1:length(MS1)
    epsi1=[epsi1;MS1{f}.epsilon(:,1)];
    epsi2=[epsi2;MS1{f}.epsilon(:,2)];
    Pr_id=[Pr_id; f+0*MS1{f}.epsilon(:,1)];
    scan_id=[scan_id 1:MS1{f}.nbscan];
end

[sepsi1,I1]=sort(epsi1);
[sepsi2,I2]=sort(epsi2);

%% get the lowest 10 scans
for i=1:10
    Pid=Pr_id(I1(i));
    sc_id=scan_id(I1(i));
    indscan=MS1{Pid}.indscan{sc_id};
    T=interp1(CTDProfiles.datadown{Pid}.ctdtime,CTDProfiles.datadown{Pid}.T,EpsiProfiles.datadown{Pid}.epsitime(indscan));
    P=interp1(CTDProfiles.datadown{Pid}.ctdtime,CTDProfiles.datadown{Pid}.P,EpsiProfiles.datadown{Pid}.epsitime(indscan));
    S=interp1(CTDProfiles.datadown{Pid}.ctdtime,CTDProfiles.datadown{Pid}.S,EpsiProfiles.datadown{Pid}.epsitime(indscan));
    scans{i}.P=P;
    scans{i}.T=T;
    scans{i}.S=S;
    scans{i}.t1=EpsiProfiles.datadown{Pid}.t1(indscan);
    scans{i}.t2=EpsiProfiles.datadown{Pid}.t2(indscan);
    scans{i}.s1=EpsiProfiles.datadown{Pid}.s1(indscan);
    scans{i}.s2=EpsiProfiles.datadown{Pid}.s2(indscan);
    scans{i}.a1=EpsiProfiles.datadown{Pid}.a1(indscan);
    scans{i}.a2=EpsiProfiles.datadown{Pid}.a2(indscan);
    scans{i}.a3=EpsiProfiles.datadown{Pid}.a3(indscan);
    scans{i}.w=MS1{Pid}.w(sc_id);
end

