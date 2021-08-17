% function process_streaming()

% when you are dome checking type:
% ctrl-C 
% Important you need to close the serial port
% delete(instrfind)


Meta_Data.epsi.s1.ADCconf='Bipolar';
Meta_Data.epsi.s2.ADCconf='Bipolar';
Meta_Data.epsi.t1.ADCconf='Unipolar';
Meta_Data.epsi.t2.ADCconf='Unipolar';
Meta_Data.epsi.a1.ADCconf='Unipolar';
Meta_Data.epsi.a2.ADCconf='Unipolar';
Meta_Data.epsi.a3.ADCconf='Unipolar';
Meta_Data.epsi.s1.full_range=2.5;
Meta_Data.epsi.s2.full_range=2.5;
Meta_Data.epsi.t1.full_range=2.5;
Meta_Data.epsi.t2.full_range=2.5;
Meta_Data.epsi.a1.full_range=1.8;
Meta_Data.epsi.a2.full_range=1.8;
Meta_Data.epsi.a3.full_range=1.8;


%% plot stuff
% function process_streaming_frompython()
addpath(genpath(fullfile('../','toolbox')));

fig=figure(1);
set(1,'units','inch','position',[0,0,35,15]);
ax(1)=subplot('Position',[.1 .89 .8 .06]);
ax(2)=subplot('Position',[.1 .81 .8 .06]);
ax(3)=subplot('Position',[.1 .73 .8 .06]);
ax(4)=subplot('Position',[.1 .65 .8 .06]);
ax(5)=subplot('Position',[.1 .24 .8 .38]);
ax(6)=subplot('Position',[.1 .08 .8 .15]);


%Ylim of the plots
alimm=-1.1;alimp=1.1;
slimm=1;slimp=2;
tlimm=1.2;tlimp=3.3;
splimm=9e-17;splimp=1e-7;

efe.namechannels={'t1','t2', ...
    's1','s2', ...
    'a1','a2','a3'};
efe.nb_channel=length(efe.namechannels);
efe.FS=325;
efe.nb_sample=160;
efe.nbblock_diag=2;
efe.dof=2;
efe.nb_block=efe.nbblock_diag*efe.dof;

% accelerometer Voltage into Accelereation units (in g).
full_range = 2.5;
full_range_accell = 1.8;
bit_counts = 24;
gain = 1;
acc_offset = 1.65;
acc_factor = 0.66;


% prep plot
cmap=colormap(parula(8));
ylabel(ax(1),'g','FontSize',20)
ylabel(ax(2),'V','FontSize',20)
ylabel(ax(3),'V','FontSize',20)
ylabel(ax(4),'sample','FontSize',20)
ylabel(ax(5),'V^2/Hz','FontSize',20)
for a=1:4
    ax(a).XTickLabel='';
    ax(a).FontSize=20;
end
ax(4).FontSize=20;
xlabel(ax(4),'(seconds)','fontsize',20)



timeaxis=linspace(0,efe.nb_block*efe.nb_sample./efe.FS,efe.nb_block*efe.nb_sample);

% spectrum stuff
% sample rate channels
efe.tscan=efe.nbblock_diag./2;
% number of samples per scan (1s) in channels
efe.df        = 1/efe.tscan;
efe.f=(efe.df:efe.df:efe.FS/2)'; % frequency vector for spectra
data=nan(efe.nb_channel,efe.dof,efe.nbblock_diag*efe.nb_sample);

Fn    = .5*efe.FS;  % Nyquist frequency
FR    = 2.5;    % Full range in Volts
def_noise=@(x)((FR/2^x)^2 /Fn);
Accelnoise=45e-6^2+0*efe.f;
set(ax(5),'fontsize',30)
ylabel(ax(5),'V^2 / Hz','fontsize',30)
xlabel(ax(5),'Hz','fontsize',30)
% title(ax(1),'coucou','fontsize',25)
grid(ax(5),'on')
% bit noise
n20=loglog(ax(5),efe.f,efe.f*0+def_noise(20),'--','Color',[.5 .5 .5],'linewidth',2);
hold(ax(5),'on')
n24=loglog(ax(5),efe.f,efe.f*0+def_noise(24),'--','Color',[.1 .1 .1],'linewidth',2);
n16=loglog(ax(5),efe.f,efe.f*0+def_noise(16),'.-','Color',[.3 .3 .3],'linewidth',2);
An=loglog(ax(5),efe.f,Accelnoise,'--','Color',[.1 .1 .1],'linewidth',2);
hold(ax(5),'off')


% now we want to open epsi_data.bin read and process the las 5 seconds
% we will plot time series and spectra with 5 seconds length.
count=0;
test=1;

% set(1,'WindowKeyPressFcn','disp(''hi''); test=0;');

fid=fopen('../STREAMING/epsi_data','r');
isread=0;
isclock=0;
while test
    if (isread==0 && isclock==0)
        tic
        isclock=1;
    end
    while isread==0
        fseek(fid,0,1);
        if ftell(fid)>5000*efe.nb_block
            isread=1;
            last_read0=ftell(fid);
            toc
        end
    end
    % 5000 is arbritrarily chosen.
    fseek(fid,-5000*efe.nb_block,1);
    str = fread(fid,'*char')';
    last_read = ftell(fid);
    if last_read>last_read0
        last_read0=last_read;
        efe.epsi=mod_som_read_epsi_raw_str(str,Meta_Data);
        
        efe.epsi =structfun(@(x) reshape(x',[],1),efe.epsi,'un',0);
        
        for c=1:efe.nb_channel
            wh_channel=efe.namechannels{c};
            tempo=reshape(efe.epsi.(wh_channel)(1:efe.dof*efe.nbblock_diag*160).',[efe.nbblock_diag*efe.nb_sample efe.dof]).';
            data(c,:,:)=tempo;
            %         [P11(c,:),f1]=pwelch(efe.epsi.(wh_channel),[],[],[],efe.FS);
        end
        % compute spectra
        
        % compute spectra
        f=efe.f;
        % compute spectra
        [f1,~,P11,Co12]=get_profile_spectrum(data,f);
        indf1=find(f1>=0);
        indf1=indf1(1:end-1);
        f1=f1(indf1);
        P11= 2*P11(:,:,indf1);
        Co12=Co12(:,:,:,indf1);
        P11_0= P11;
        P11 = smoothdata(P11,3,'gaussian',5);
        
        % set coheence correction
        inda1=5;inda2=6;inda3=7;
        indt1=1;indt2=2;
        inds1=3;inds2=4;
        a1f=nanmean(squeeze(P11_0(inda1,:,:)));
        a2f=nanmean(squeeze(P11_0(inda2,:,:)));
        a3f=nanmean(squeeze(P11_0(inda3,:,:)));
        
        s1=squeeze(P11_0(inds1,:,:));
        s1f=nanmean(s1);
        Cos1a3=squeeze(Co12(inds1,inda3-1,:,:));
        Cos1a2=squeeze(Co12(inds1,inda2-1,:,:));
        Cos1a1=squeeze(Co12(inds1,inda1-1,:,:));
        
        Cos1a3=abs(nanmean(Cos1a3)).^2./s1f./a3f;
        Cos1a2=abs(nanmean(Cos1a2)).^2./s1f./a2f;
        Cos1a1=abs(nanmean(Cos1a1)).^2./s1f./a1f;
        %         Cos1tot=max(cat(3,Cos1a1,Cos1a2,Cos1a3),[],3);
        Cos1tot=smoothdata(Cos1a3,'gaussian',3);
        ib=find(Cos1tot>1);  Cos1tot(ib)=1;  %TEMPTEMPTEMP fix
        
        s2=squeeze(P11_0(inds2,:,:));
        s2f=nanmean(s2);
        Cos2a3=squeeze(Co12(inds2,inda3-1,:,:));
        Cos2a2=squeeze(Co12(inds2,inda2-1,:,:));
        Cos2a1=squeeze(Co12(inds2,inda1-1,:,:));
        Cos2a3=abs(nanmean(Cos2a3)).^2./s2f./a3f;
        Cos2a2=abs(nanmean(Cos2a2)).^2./s2f./a2f;
        Cos2a1=abs(nanmean(Cos2a1)).^2./s2f./a1f;
        Cos2tot=smoothdata(Cos2a3,'gaussian',3);
        ib=find(Cos2tot>1);  Cos2tot(ib)=1;  %TEMPTEMPTEMP fix
        s1=squeeze(P11_0(inds1,:,:)); s1=s1.*(1-Cos1tot);
        s2=squeeze(P11_0(inds2,:,:)); s2=s2.*(1-Cos2tot);
        s1 = smoothdata(s1,2,'gaussian',5);
        s2 = smoothdata(s2,2,'gaussian',5);
        t1=squeeze(P11_0(indt1,:,:));
        t2=squeeze(P11_0(indt2,:,:));
        t1 = smoothdata(t1,2,'gaussian',5);
        t2 = smoothdata(t2,2,'gaussian',5);
        
        MS.Pf_0=P11;
        MS.f=f1;
        
        % plot time series
        plot(ax(1),timeaxis,detrend(efe.epsi.a1(1:efe.nb_block*efe.nb_sample)),'Color',cmap(1,:),'linewidth',2)
        hold(ax(1),'on')
        plot(ax(1),timeaxis,detrend(efe.epsi.a2(1:efe.nb_block*efe.nb_sample)),'Color',cmap(2,:),'linewidth',2)
        plot(ax(1),timeaxis,detrend(efe.epsi.a3(1:efe.nb_block*efe.nb_sample)),'Color',cmap(3,:),'linewidth',2)
        hold(ax(1),'off')
        
        legend(ax(1),{sprintf('a1=%1.2fV',nanmean(efe.epsi.a1(1:efe.nb_block*efe.nb_sample))), ...
            sprintf('a2=%1.2fV',nanmean(efe.epsi.a2(1:efe.nb_block*efe.nb_sample))), ...
            sprintf('a3=%1.2fV',nanmean(efe.epsi.a3(1:efe.nb_block*efe.nb_sample)))})
        
        plot(ax(2),timeaxis,detrend(efe.epsi.t1(1:efe.nb_block*efe.nb_sample)),'Color',cmap(4,:),'linewidth',2)
        hold(ax(2),'on')
        plot(ax(2),timeaxis,detrend(efe.epsi.t2(1:efe.nb_block*efe.nb_sample)),'Color',cmap(5,:),'linewidth',2)
        hold(ax(2),'off')
        legend(ax(2),{sprintf('t1=%1.2fV',nanmean(efe.epsi.t1(1:efe.nb_block*efe.nb_sample))), ...
            sprintf('t2=%1.2fV',nanmean(efe.epsi.t2(1:efe.nb_block*efe.nb_sample)))})
        
        
        
        plot(ax(3),timeaxis,detrend(efe.epsi.s1(1:efe.nb_block*efe.nb_sample)),'Color',cmap(6,:),'linewidth',2)
        hold(ax(3),'on')
        plot(ax(3),timeaxis,detrend(efe.epsi.s2(1:efe.nb_block*efe.nb_sample)),'Color',cmap(7,:),'linewidth',2)
        hold(ax(3),'off')
        legend(ax(3),{sprintf('s1=%1.2fV',nanmean(efe.epsi.s1(1:efe.nb_block*efe.nb_sample))), ...
            sprintf('s2=%1.2fV',nanmean(efe.epsi.s2(1:efe.nb_block*efe.nb_sample)))})
        
        %         plot(ax(4),timeaxis(1:end-1),diff(efe.epsi.c),'Color',cmap(8,:),'linewidth',2)
        
        %         legend(ax(1),{'a1','a2','a3'})
        %         legend(ax(2),{'t1','t2'})
        %         legend(ax(3),{'s1','s2'})
        legend(ax(4),{'diff ramp'})
        
        %%
        ind_5Hz35Hz=find(MS.f>5 & MS.f<35);
        a3_5Hz35Hz=squeeze(nanmean(nansum(MS.Pf_0(inda3,:,ind_5Hz35Hz),3)))*MS.f(2);
        s1_5Hz35Hz=squeeze(nanmean(nansum(s1(:,ind_5Hz35Hz),2)))*MS.f(2);
        s2_5Hz35Hz=squeeze(nanmean(nansum(s2(:,ind_5Hz35Hz),2)))*MS.f(2);
        t1_5Hz35Hz=squeeze(nanmean(nansum(t1(:,ind_5Hz35Hz),2)))*MS.f(2);
        t2_5Hz35Hz=squeeze(nanmean(nansum(t2(:,ind_5Hz35Hz),2)))*MS.f(2);
        %%
        % plot spectra
        hold(ax(5),'on')
        l0=loglog(ax(5),f1,squeeze(nanmean(P11(indt1,:,:),2)),'Color',cmap(4,:),'linewidth',2);
        l1=loglog(ax(5),f1,squeeze(nanmean(P11(indt2,:,:),2)),'Color',cmap(5,:),'linewidth',2);
        l21=loglog(ax(5),f1,squeeze(nanmean(P11(inds1,:,:),2)),'--','Color',cmap(6,:),'linewidth',2);
        l31=loglog(ax(5),f1,squeeze(nanmean(P11(inds2,:,:),2)),'--','Color',cmap(7,:),'linewidth',2);
        l2=loglog(ax(5),f1,squeeze(nanmean(s1)),'Color',cmap(6,:),'linewidth',2);
        l3=loglog(ax(5),f1,squeeze(nanmean(s2)),'Color',cmap(7,:),'linewidth',2);
        %         l4=loglog(ax(5),f1,squeeze(nanmean(P11(5,:,:),2)),'Color',cmap(1,:));
        l4=loglog(ax(5),f1,squeeze(nanmean(P11(inda1,:,:),2)),'Color',cmap(1,:),'linewidth',2);
        l5=loglog(ax(5),f1,squeeze(nanmean(P11(inda2,:,:),2)),'Color',cmap(2,:),'linewidth',2);
        %         l6=loglog(ax(5),f1,squeeze(nanmean(P11(8,:,:),2)),'Color',cmap(3,:),'linewidth',2);
        l6=loglog(ax(5),f1,a3f,'Color',cmap(3,:),'linewidth',2);
        set(ax(5),'Xscale','log','Yscale','log')
        legend([n24 n20 n16 An l0 l1 l2 l3 l4 l5 l6],{'24 bit','20 bit','16 bit','Accel noise','t1','t2','s1','s2','a1','a2','a3'},'location','SouthWest')
        F=fill(ax(5),[f1(ind_5Hz35Hz([1 end])) f1(ind_5Hz35Hz([end 1]))],[1e-16 1e-16 1e-4 1e-4],'k');
        F.FaceAlpha=.1;
        T1=text(1,1e-6*.6,sprintf('s1=%2.2eV',sqrt(s1_5Hz35Hz)),'fontsize',25,'Parent',ax(5));
        T2=text(1,1e-6*.1,sprintf('s2=%2.2eV',sqrt(s2_5Hz35Hz)),'fontsize',25','Parent',ax(5));
        T3=text(1,1e-6*.006,sprintf('t1=%2.2eV',sqrt(t1_5Hz35Hz)),'fontsize',25,'Parent',ax(5));
        T4=text(1,1e-6*.001,sprintf('t2=%2.2eV',sqrt(t2_5Hz35Hz)),'fontsize',25','Parent',ax(5));
        
        hold(ax(5),'off')
        
        semilogx(ax(6),f1,Cos1tot,'Color',cmap(6,:),'linewidth',2)
        hold(ax(6),'on')
        semilogx(ax(6),f1,Cos2tot,'Color',cmap(7,:),'linewidth',2)
        hold(ax(6),'off')
        
        ax(5).XLim=[efe.df f(end)];
        ax(6).XLim=[efe.df f(end)];
        
        ax(1).YLim=.001*[-1 1];
        ax(2).YLim=.00001*[-1 1];
        ax(3).YLim=.00005*[-1 1];
        ax(4).YLim=[0 2];
        ax(6).YLim=[0 1];
        ax(5).YLim=[splimm splimp];
        ax(1).XLim=[0 4];
        ax(2).XLim=[0 4];
        ax(3).XLim=[0 4];
        ax(4).XLim=[0 4];
        for p=1:6
            ax(p).FontSize=20;
            if p<4
                ax(p).XTickLabel='';
            end
            if p==5
                ax(p).XTickLabel='';
            end
        end
        
        ylabel(ax(1),'V','FontSize',20)
        ylabel(ax(2),'V','FontSize',20)
        ylabel(ax(3),'V','FontSize',20)
        ylabel(ax(4),'sample','FontSize',20)
        ylabel(ax(5),'V^2/Hz','FontSize',20)
        xlabel(ax(6),'Hz','FontSize',20)
        
        
        
        pause(.001)
        delete(l0);
        delete(l1);
        delete(l2);
        delete(l3);
        delete(l21);
        delete(l31);
        delete(l4);
        delete(l5);
        delete(l6);
        delete(F);
        delete(T1);
        delete(T2);
        delete(T3);
        delete(T4);
    else
        fprintf('last read %i\n',last_read)
        fprintf('last read0 %i\n',last_read0)
        isread=0;
        isclock=0;
    end
end

function [k,P1,P11,Co12]=get_profile_spectrum(data,k)
%
%  input: data
% . data : epsi data
% . k:  frequency array
%  Created by Arnaud Le Boyer on 7/28/18.


switch length(size(data))
    case 3 % reshape into 2D matrice for fft and then reshape
        [nb_sensor,nb_scan,Lscan]=size(data);
        data=reshape(data,[nb_sensor* nb_scan Lscan]);
        Lax1=nb_sensor* nb_scan;
        size_data=3;
    case 2
        [Lax1,Lscan]=size(data);
        size_data=2;
    otherwise
        warning('no valid size for data : get power spectrum')
end

dk=k(1);
window = ones(Lax1,1)*hamming(Lscan).';
wc2=1/mean(window(1,:).^2);            % window correction factor
% datap  = window.*(data- mean(data,2)* ones(1,Lscan));
datap=window.* detrend(data.').';
P1  = fft(datap,[],2);
P11 = conj(P1).*P1./Lscan^2/dk*wc2;

if size_data==3
    P1=reshape(P1,[nb_sensor,nb_scan,Lscan]);
    P11=reshape(P11,[nb_sensor,nb_scan,Lscan]);
    P12=zeros(nb_sensor,nb_sensor-1,nb_scan,Lscan);
    Co12=zeros(nb_sensor,nb_sensor-1,nb_scan,Lscan);
    ind_nbsensor=1:nb_sensor;
    for j=1:nb_sensor
        tempo=shiftdim(repmat(squeeze(P1(j,:,:)),[1,1,nb_sensor-1]),2);
        P12(j,:,:,:)  = conj(tempo).*P1(ind_nbsensor~=j,:,:)./Lscan^2/dk*wc2;
        Co12(j,:,:,:) = 2*squeeze(P12(j,:,:,:));
    end
end

if rem(Lscan,2)==0
    k=-Lscan/2*dk:dk:Lscan/2*dk-dk;
else
    kp=dk:dk:dk*floor(Lscan/2);
    k=[fliplr(-kp) 0 kp];
end
k=fftshift(k);


end

