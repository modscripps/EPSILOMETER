function EPSI_matchCTDtime(Meta_Data,ctdfile)

%  When epsi on WW. CTD and epsi are not on the time. 
% . We are looking to reconcile the indexes
%
%  input: Meta_Data
%  created with Meta_Data=create_Meta_Data(file). Meta_Data contain the
%  path to calibration file and EPSI configuration needed to process the
%  epsi data
%
%  Created by Arnaud Le Boyer on 7/28/18.

ctdfile='offset_raw.mat';
ctdname=ctdfile(1:end-4);


CTDpath=Meta_Data.CTDpath;
Epsipath=Meta_Data.Epsipath;

CTD=load(fullfile(CTDpath,ctdfile));
t1=[];
load(fullfile(Epsipath,['epsi_' Meta_Data.deployement]),...
                        'epsitime','t1','t2','s1','s2','c','a1','a2','a3')

%%
% datestrctd1='06-Jun-2018 06:00:23';
% indctd=find(CTD.ctdtime==ctd1.Position(1));
indctd=208777;  % MISOBoB 2018
epsi1.Position(1)=3144178;% MISOBoB 2018

fsampling=325.5;

P=CTD.(ctdname).P;
T=CTD.(ctdname).T;
S=CTD.(ctdname).S;
C=CTD.(ctdname).C;
sig=CTD.(ctdname).sig;
ctdtime=CTD.(ctdname).time;

epsitime1=ctdtime(indctd)+ (0:length(t1)-epsi1.Position(1))./fsampling/86400;
t1 = t1(epsi1.Position(1):end);
t2 = t2(epsi1.Position(1):end);
s1 = s1(epsi1.Position(1):end);
s2 = s2(epsi1.Position(1):end);
c  = c(epsi1.Position(1):end);
a1 = a1(epsi1.Position(1):end);
a2 = a2(epsi1.Position(1):end);
a3 = a3(epsi1.Position(1):end);

close all
ax(1)=subplot(211);
plot(epsitime1,t1)
ax(2)=subplot(212);
plot(ctdtime,T)
xlim(epsitime1([1 end]))

linkaxes(ax,'x')


epsitime=epsitime1;


save(fullfile(CTDpath,['ctd_' Meta_Data.deployement]),'ctdtime','P','T','C','S','sig')
save(fullfile(Epsipath,['epsi_' Meta_Data.deployement]),'epsitime','t1','t2','s1','s2','c','a1','a2','a3')

