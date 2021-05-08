% plot time series and fft (Pf_0) from PISTON 2019
%process_dir='~/ARNAUD/SCRIPPS/EPSILOMETER/';
addpath(process_dir)
EPSI_matlab_path(process_dir); % EPSI_matlab_path is under ~/ARNAUD/SCRIPPS/EPSILOMETER/
rem_nan=@(x) (fillmissing(x,'linear'));

load('/Volumes/Ahua/data_archive/WaveChasers-DataArchive/PISTON/Cruises/sr1914/share/data/epsi/PC2/d12/L1/Turbulence_Profiles0.mat')
load('/Volumes/Ahua/data_archive/WaveChasers-DataArchive/PISTON/Cruises/sr1914/share/data/epsi/PC2/d12/L1/Profiles_d10.mat')
load('/Volumes/Ahua/data_archive/WaveChasers-DataArchive/PISTON/Cruises/sr1914/share/data/epsi/PC2/d12/L1/Meta_Data.mat')
% load('/Volumes/DataDrive/NISKINE/NISKINE/epsifish1bis/d4/L1/Meta_Data.mat')
% load('/Volumes/DataDrive/NISKINE/NISKINE/epsifish1bis/d4/L1/Turbulence_Profiles0.mat')
% load('/Volumes/DataDrive/NISKINE/NISKINE/epsifish1bis/d4/L1/Profiles_d4.mat')

id_profile=1;



% All_channels={'t1','t2','s1','s2','c','a1','a2','a3'};
All_channels=Meta_Data.PROCESS.channels;
nb_channels=length(All_channels);

nbscan  = MS{id_profile}.nbscan;
Lscan   = length(MS{id_profile}.indscan{1});

try
CTD_Profiles=CTDProfiles.datadown;
EPSI_Profiles=EpsiProfiles.datadown;
catch
CTD_Profiles=CTDProfile.datadown;
EPSI_Profiles=EpsiProfile.datadown;
end


Profile=EPSI_Profiles{id_profile};
Profile.P=interp1(rem_nan(CTD_Profiles{id_profile}.ctdtime),CTD_Profiles{id_profile}.P,Profile.epsitime);
Profile.P=filloutliers(Profile.P,'center','movmedian',1000);
Profile.T=interp1(rem_nan(CTD_Profiles{id_profile}.ctdtime),CTD_Profiles{id_profile}.T,Profile.epsitime);
Profile.S=interp1(rem_nan(CTD_Profiles{id_profile}.ctdtime),CTD_Profiles{id_profile}.S,Profile.epsitime);



% split time serie to make data matrice
dataarray=zeros(nb_channels,nbscan,Lscan);

for c=1:length(All_channels)
    wh_channels=All_channels{c};
    ind=find(cellfun(@(x) strcmp(x,wh_channels),All_channels));
    switch wh_channels
        case {'t1','t2','s1','s2'}
            data.(wh_channels) = cell2mat(cellfun(@(x) filloutliers(Profile.(wh_channels)(x),'center','movmedian',ceil(MS{id_profile}.f(end))),MS{id_profile}.indscan,'un',0)).';
            dataclean.(wh_channels) = cell2mat(cellfun(@(x) Profile.(wh_channels)(x),MS{id_profile}.indscan,'un',0)).';
            dataarray(ind,:,:) = dataclean.(wh_channels); 
        case {'a1','a2','a3'}
            data.(wh_channels) = cell2mat(cellfun(@(x) filloutliers(Profile.(wh_channels)(x),'center','movmedian',ceil(MS{id_profile}.f(end))),MS{id_profile}.indscan,'un',0)).';
            dataclean.(wh_channels) = cell2mat(cellfun(@(x) Profile.(wh_channels)(x),MS{id_profile}.indscan,'un',0)).';
            dataarray(ind,:,:) = dataclean.(wh_channels); 
        case {'c'}
            data.(wh_channels) = cell2mat(cellfun(@(x) Profile.(wh_channels)(x),MS{id_profile}.indscan,'un',0)).';
            dataarray(ind,:,:) = data.(wh_channels); 
    end
end

% compute spectra
tscan     =  3;
FS        = 320;
df        = 1/tscan;
f=(df:df:FS/2)'; % frequency vector for spectra
[f1,~,P11,~]=get_profile_spectrum(dataarray,f);
indf1=find(f1>=0);
indf1=indf1(1:end-1);
f1=f1(indf1);
P11= 2*P11(:,:,indf1);



%%
tscan=length(MS{id_profile}.indscan{1})/320;
df=1./tscan;
timeaxis=linspace(0,tscan,tscan*320);


close all


figure('units','inch','position',[0,0,35,15]);
ax(1)=subplot('Position',[.3 .83 .6 .1]);
ax(2)=subplot('Position',[.3 .72 .6 .1]);
ax(3)=subplot('Position',[.3 .61 .6 .1]);
ax(4)=subplot('Position',[.3 .5 .6 .1]);
ax(5)=subplot('Position',[.3 .05 .6 .4]);
ax(6)=subplot('Position',[.1 .1 .15 .7]);

semilogx(ax(6),MS{id_profile}.epsilon,MS{id_profile}.pr)
cmap=colormap(parula(8));

%Ylim of the plots
alimm=-1.1;alimp=1.1;
slimm=1.2;slimp=1.3;
tlimm=1.2;tlimp=3.3;
splimm=9e-17;splimp=1e-4;

% bit noise level
FS=320;
Fn    = .5*FS;  % Nyquist frequency
FR    = 2.5;    % Full range in Volts
def_noise=@(x)((FR/2^x)^2 /Fn);
Accelnoise=45e-6^2+0*MS{id_profile}.f;
set(ax(5),'fontsize',30)
ylabel(ax(5),'V^2 / Hz','fontsize',30)
xlabel(ax(5),'Hz','fontsize',30)
title(ax(1),'coucou','fontsize',25)
grid(ax(5),'on')
% bit noise
n20=loglog(ax(5),MS{id_profile}.f,MS{id_profile}.f*0+def_noise(20),'--','Color',[.5 .5 .5],'linewidth',2);
hold(ax(5),'on')
n24=loglog(ax(5),MS{id_profile}.f,MS{id_profile}.f*0+def_noise(24),'--','Color',[.1 .1 .1],'linewidth',2);
n16=loglog(ax(5),MS{id_profile}.f,MS{id_profile}.f*0+def_noise(16),'.-','Color',[.3 .3 .3],'linewidth',2);
An=loglog(ax(5),MS{id_profile}.f,Accelnoise,'--','Color',[.1 .1 .1],'linewidth',2);
% loglog(ax(5),MS{id_profile}.f,nanmedian(squeeze(MS{id_profile}.Pf_0(1,:,:))),'g');
loglog(ax(5),MS{id_profile}.f,nanmedian(squeeze(MS{id_profile}.Pf_0(3,:,:))),'m');
hold(ax(5),'off')

v = VideoWriter(sprintf('%s_rawcast%i%s','NISKINE_epsi1_d4',id_profile,'.avi'));
v.FrameRate=5;
open(v)


for idscan=1:MS{id_profile}.nbscan

    semilogx(ax(6),MS{id_profile}.epsilon(idscan,:),MS{id_profile}.pr(idscan),'kp')

% plot time series
plot(ax(1),timeaxis,detrend(data.a1(idscan,:)),'Color',cmap(1,:))
hold(ax(1),'on')
plot(ax(1),timeaxis,detrend(data.a2(idscan,:)),'Color',cmap(2,:))
plot(ax(1),timeaxis,detrend(data.a3(idscan,:)),'Color',cmap(3,:))
hold(ax(1),'off')
title(ax(1),sprintf('scan %i over %i',idscan,nbscan))

plot(ax(2),timeaxis,data.t1(idscan,:),'Color',cmap(4,:))
hold(ax(2),'on')
plot(ax(2),timeaxis,dataclean.t1(idscan,:),'--','Color',cmap(4,:))
plot(ax(2),timeaxis,data.t2(idscan,:),'Color',cmap(5,:))
plot(ax(2),timeaxis,dataclean.t2(idscan,:),'--','Color',cmap(5,:))
hold(ax(2),'off')

plot(ax(3),timeaxis,data.s1(idscan,:),'Color',cmap(6,:))
hold(ax(3),'on')
plot(ax(3),timeaxis,dataclean.s1(idscan,:),'--','Color',cmap(6,:))
plot(ax(3),timeaxis,data.s2(idscan,:),'Color',cmap(7,:))
plot(ax(3),timeaxis,dataclean.s2(idscan,:),'--','Color',cmap(7,:))
hold(ax(3),'off')

plot(ax(4),timeaxis(1:end-1),diff(data.c(idscan,:)),'Color',cmap(8,:))

legend(ax(1),{'detrend(a1)','detrend(a2)','detrend(a3)'})
legend(ax(2),{'t1','t2'})
legend(ax(3),{'s1','s2'})
legend(ax(4),{'diff ramp'})


% plot sepctra

hold(ax(5),'on')
% l0=loglog(ax(5),MS{id_profile}.f,squeeze(MS{id_profile}.Pf_0(1,idscan,:)),'Color',cmap(4,:));
% l1=loglog(ax(5),MS{id_profile}.f,squeeze(MS{id_profile}.Pf_0(2,idscan,:)),'Color',cmap(5,:));
l2=loglog(ax(5),MS{id_profile}.f,squeeze(MS{id_profile}.Pf_0(3,idscan,:)),'Color',cmap(6,:));
l3=loglog(ax(5),MS{id_profile}.f,squeeze(MS{id_profile}.Pf_0(4,idscan,:)),'Color',cmap(7,:));
l22=loglog(ax(5),MS{id_profile}.f,squeeze(P11(3,idscan,:)),'+-','Color',cmap(6,:));
l33=loglog(ax(5),MS{id_profile}.f,squeeze(P11(4,idscan,:)),'+-','Color',cmap(7,:));
%         l4=loglog(ax(5),MS{id_profile}.f,squeeze(nanmean(MS{id_profile}.Pf_0(5,:,:),2)),'Color',cmap(1,:));
% l4=loglog(ax(5),MS{id_profile}.f,squeeze(MS{id_profile}.Pf_0(6,idscan,:)),'Color',cmap(1,:));
% l5=loglog(ax(5),MS{id_profile}.f,squeeze(MS{id_profile}.Pf_0(7,idscan,:)),'Color',cmap(2,:));
% l6=loglog(ax(5),MS{id_profile}.f,squeeze(MS{id_profile}.Pf_0(8,idscan,:)),'Color',cmap(3,:));
set(ax(5),'Xscale','log','Yscale','log')
% legend([n24 n20 n16 An l0 l1 l2 l3 l4 l5 l6],{'24 bit','20 bit','16 bit','Accel noise','t1','t2','s1','s2','a1(g^2/Hz)','a2(g^2/Hz)','a3(g^2/Hz)'},'location','SouthWest')
legend([n24 n20 n16 An l2 l3],{'24 bit','20 bit','16 bit','Accel noise','s1','s2'},'location','SouthWest')
hold(ax(5),'off')
ax(5).XLim=[df ceil(MS{id_profile}.f(end))];

% ax(1).YLim=[alimm alimp];
% ax(2).YLim=[slimm slimp];
% ax(3).YLim=[tlimm tlimp];
% ax(4).YLim=[0 2];
ax(5).YLim=[splimm splimp];
ax(1).XLim=[0 3];
ax(2).XLim=[0 3];
ax(3).XLim=[0 3];
ax(4).XLim=[0 3];
for p=1:5
    ax(p).FontSize=20;
end

ylabel(ax(1),'g','FontSize',20)
ylabel(ax(2),'V','FontSize',20)
ylabel(ax(3),'V','FontSize',20)
ylabel(ax(4),'sample','FontSize',20)
ylabel(ax(5),'V^2/Hz','FontSize',20)

% 
frame=getframe(gcf);
writeVideo(v,frame)
pause(.001)
 
% delete(l0);
% delete(l1);
delete(l2);
delete(l3);
delete(l22);
delete(l33);
% delete(l4);
% delete(l5);
% delete(l6);


end

close(v)












