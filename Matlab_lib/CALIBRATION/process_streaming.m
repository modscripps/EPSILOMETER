% function process_streaming()

% when you are dome checking type:
% ctrl-C 
% Important you need to close the serial port
% delete(instrfind)



% open serial port 
listserial=seriallist;
for l=1:numel(listserial)
    fprintf('%i=%s\n',l,listserial{l})
end
id_port=input('Choose your port\nPort:');

 ser=serial(listserial{id_port},'BaudRate',115200);
%ser=serial(listserial{6},'BaudRate',115200);
fopen(ser);

try 
    %% plot stuff
    addpath(genpath(fullfile('..','toolbox')));
    
    figure('units','inch','position',[0,0,35,15]);
    ax(1)=subplot('Position',[.1 .83 .8 .1]);
    ax(2)=subplot('Position',[.1 .72 .8 .1]);
    ax(3)=subplot('Position',[.1 .61 .8 .1]);
    ax(4)=subplot('Position',[.1 .5 .8 .1]);
    ax(5)=subplot('Position',[.1 .05 .8 .4]);
    
    %Ylim of the plots
    alimm=-1.1;alimp=1.1;
    slimm=0;slimp=2.5;
    tlimm=0;tlimp=2.5;
    splimm=9e-17;splimp=1e-6;
    
    % define some important and fixed variable
    epsi.bytes_per_channel=3; % ADC is 3 bytes
    epsi.nbsamples=160;% number of epsi blocks is 160 ~ 0.5 seconds
    epsi.nbblock_diag=10;% 10*.5sec blocks = 5 seconds
    epsi.name_length=5; % 5 bytes EPSI
    epsi.finishblock=2; % 2 bytes \r\n
    namechannels={'t1','t2', ...
        's1','s2', ...
        'c', ...
        'a1','a2','a3'};
    countconversion={'Unipolar','Unipolar', ...
        'Unipolar','Unipolar', ...
        'Unipolar', ...
        'Unipolar','Unipolar','Unipolar'};
    
    % accelerometer Voltage into Accelereation units (in g).
    full_range = 2.5;
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
    
    
    
    length_diag=(epsi.nbblock_diag)*epsi.nbsamples; % 1760 sample
    timeaxis=linspace(0,length_diag/325,length_diag);
    
    % spectrum stuff
    % sample rate channels
    nb_segment=5;
    FS        = 325;
    tscan=5.5;
    % number of samples per scan (1s) in channels
    df        = 1/tscan;
    f=(df:df:FS/2)'; % frequency vector for spectra
    data=nan(length(namechannels),nb_segment,length_diag);
    Fn    = .5*FS;  % Nyquist frequency
    FR    = 2.5;    % Full range in Volts
    def_noise=@(x)((FR/2^x)^2 /Fn);
    Accelnoise=45e-6^2+0*f;
    set(ax(5),'fontsize',30)
    ylabel(ax(5),'V^2 / Hz','fontsize',30)
    xlabel(ax(5),'Hz','fontsize',30)
    title(ax(1),'coucou','fontsize',25)
    grid(ax(5),'on')
    % bit noise
    n20=loglog(ax(5),f,f*0+def_noise(20),'--','Color',[.5 .5 .5],'linewidth',2);
    hold(ax(5),'on')
    n24=loglog(ax(5),f,f*0+def_noise(24),'--','Color',[.1 .1 .1],'linewidth',2);
    n16=loglog(ax(5),f,f*0+def_noise(16),'.-','Color',[.3 .3 .3],'linewidth',2);
    An=loglog(ax(5),f,Accelnoise,'--','Color',[.1 .1 .1],'linewidth',2);
    hold(ax(5),'off')
    
    
    % now we want to open epsi_data.bin read and process the las 5 seconds
    % we will plot time series and spectra with 5 seconds length.
    count=0;
    laps=0;
    begin=0;
    str=[];
    while begin==0
        bytes = fscanf(ser,'%c',512);
        str=[str bytes];
        ind_madre1 = strfind(str,'$MADRE');
        if numel(ind_madre1)==2
            begin=1;
            str1=str(ind_madre1(end):end);
            blocksize=diff(ind_madre1)-1;
            disp(str1);
        end
    end
    str = char(ones(1,10*blocksize) * '0');
    tempo_str=char(ones(1,2*blocksize) * '0');
    nbbytes=numel(str1);
    tempo_str(1:numel(str1))=str1;
    counttoto=0;
    while 1
        %     fid=fopen('../STREAMING/epsi_data.bin','r');
        %     str = fread(fid,'*char')';
        counttoto1=0;
        while nbbytes<=blocksize+1
            bytes = fscanf(ser,'%c',512)
            tempo_str(nbbytes+1:nbbytes+numel(bytes))=bytes;
            nbbytes=nbbytes+numel(bytes);
            counttoto1=counttoto1+1;
        end
        disp(counttoto1)
        
        ind_madre1 = strfind(tempo_str,'$MADRE');
        str(end-blocksize+1:end)=tempo_str(1:ind_madre1(end)-2);
        str=circshift(str,-blocksize);
        %str(end-blocksize:end)
        ind_madre = strfind(str,'$MADRE');
        ind_epsi = strfind(str,'$EPSI');
        tempo_str(1:blocksize-1)=tempo_str(ind_madre1(end):end);
        nbbytes=nbbytes-blocksize-1;
        disp(counttoto)
        %     fclose(fid);
        count=mod(count,nb_segment)+1;
        
        % define some diagnostic variable
        epsi.offset=unique(ind_epsi(1:length(ind_epsi))-ind_madre(1:length(ind_epsi)))-epsi.finishblock;
        if numel(epsi.offset)>1
            warning('issue in the Headers. They are not the same length')
        end
        
        epsi.total_length=unique(diff(ind_epsi))-epsi.finishblock-epsi.name_length-epsi.offset-1;
        if numel(epsi.total_length)>1
            warning('issue in the epsi blocks. They are not the same length')
        end
        
        epsi.nchannels=epsi.total_length/epsi.bytes_per_channel/epsi.nbsamples;
        if rem(epsi.nchannels,1)>0
            warning('issue in the epsi blocks. They are not the same length')
        end
        
        if numel(ind_epsi)==10
            %convert 3 bytes ADC samples into 24 bits counts.
            epsi.raw = int32(zeros(epsi.nbblock_diag,epsi.total_length));
            %epsi.raw = cell2mat(arrayfun(@(x) int32(str(x+epsi.name_length-1+(1:epsi.total_length))),ind_epsi(end-epsi.nbblock_diag-1:end-1),'un',0).');
            epsi.raw = cell2mat(arrayfun(@(x) int32(str(x+epsi.name_length-1+(1:epsi.total_length))),ind_epsi,'un',0).');
            epsi.raw1 = epsi.raw(:,1:epsi.bytes_per_channel:end)*256^2+ ...
                epsi.raw(:,2:epsi.bytes_per_channel:end)*256+ ...
                epsi.raw(:,3:epsi.bytes_per_channel:end);
            
            % convert count in volts
            for cha=1:epsi.nchannels
                wh_channel=namechannels{cha};
                if ~strcmp(wh_channel,'c')
                    switch countconversion{cha}
                        case 'Bipolar'
                            EPSI.epsi.(wh_channel)=full_range/gain* ...
                                (double(epsi.raw1(:,cha:epsi.nchannels:end))/2.^(bit_counts-1)-1);
                        case 'Unipolar'
                            EPSI.epsi.(wh_channel)=full_range/gain* ...
                                double(epsi.raw1(:,cha:epsi.nchannels:end))/2.^(bit_counts);
                    end
                    switch wh_channel
                        case 'a1'
                            EPSI.epsi.a1 = (EPSI.epsi.a1-acc_offset)/acc_factor;
                        case 'a2'
                            EPSI.epsi.a2 = (EPSI.epsi.a2-acc_offset)/acc_factor;
                        case 'a3'
                            EPSI.epsi.a3 = (EPSI.epsi.a3-acc_offset)/acc_factor;
                    end
                    
                else
                    EPSI.epsi.(wh_channel)=double(epsi.raw1(:,cha:epsi.nchannels:end));
                end
            end

            EPSI.epsi =structfun(@(x) reshape(x',[],1),EPSI.epsi,'un',0);
            
            % compute spectra
            data(:,count,:)=struct2array(EPSI.epsi).';
            
            % compute spectra
            [f1,~,P11,~]=get_profile_spectrum(data,f);
            indf1=find(f1>=0);
            indf1=indf1(1:end-1);
            f1=f1(indf1);
            P11= 2*P11(:,:,indf1);
            
            
            % plot time series
            plot(ax(1),timeaxis,EPSI.epsi.a1,'Color',cmap(1,:))
            hold(ax(1),'on')
            plot(ax(1),timeaxis,EPSI.epsi.a2,'Color',cmap(2,:))
            plot(ax(1),timeaxis,EPSI.epsi.a3,'Color',cmap(3,:))
            hold(ax(1),'off')
            
            plot(ax(2),timeaxis,EPSI.epsi.t1,'Color',cmap(4,:))
            hold(ax(2),'on')
            plot(ax(2),timeaxis,EPSI.epsi.t2,'Color',cmap(5,:))
            hold(ax(2),'off')
            
            plot(ax(3),timeaxis,EPSI.epsi.s1,'Color',cmap(6,:))
            hold(ax(3),'on')
            plot(ax(3),timeaxis,EPSI.epsi.s2,'Color',cmap(7,:))
            hold(ax(3),'off')
            
            plot(ax(4),timeaxis(1:end-1),diff(EPSI.epsi.c),'Color',cmap(8,:))
            
            legend(ax(1),{'a1','a2','a3'})
            legend(ax(2),{'t1','t2'})
            legend(ax(3),{'s1','s2'})
            legend(ax(4),{'diff ramp'})
            
            
            % plot sepctra
            
            hold(ax(5),'on')
            l0=loglog(ax(5),f1,squeeze(nanmean(P11(1,:,:),2)),'Color',cmap(6,:));
            l1=loglog(ax(5),f1,squeeze(nanmean(P11(2,:,:),2)),'Color',cmap(7,:));
            l2=loglog(ax(5),f1,squeeze(nanmean(P11(3,:,:),2)),'Color',cmap(4,:));
            l3=loglog(ax(5),f1,squeeze(nanmean(P11(4,:,:),2)),'Color',cmap(5,:));
            l4=loglog(ax(5),f1,squeeze(nanmean(P11(5,:,:),2)),'Color',cmap(1,:));
            l5=loglog(ax(5),f1,squeeze(nanmean(P11(6,:,:),2)),'Color',cmap(2,:));
            l6=loglog(ax(5),f1,squeeze(nanmean(P11(7,:,:),2)),'Color',cmap(3,:));
            set(ax(5),'Xscale','log','Yscale','log')
            legend([n24 n20 n16 An l0 l1 l2 l3 l4 l5 l6],{'24 bit','20 bit','16 bit','Accel noise','t1','t2','s1','s2','a1(g^2/Hz)','a2(g^2/Hz)','a3(g^2/Hz)'},'location','SouthWest')
            hold(ax(5),'off')
            ax(5).XLim=[df f(end)];
            
            ax(1).YLim=[alimm alimp];
            ax(2).YLim=[slimm slimp];
            ax(3).YLim=[tlimm tlimp];
            ax(4).YLim=[0 2];
            ax(5).YLim=[splimm splimp];
            ax(1).XLim=[0 4];
            ax(2).XLim=[0 4];
            ax(3).XLim=[0 4];
            ax(4).XLim=[0 4];
            for p=1:5
                ax(p).FontSize=20;
            end
            
            ylabel(ax(1),'g','FontSize',20)
            ylabel(ax(2),'V','FontSize',20)
            ylabel(ax(3),'V','FontSize',20)
            ylabel(ax(4),'sample','FontSize',20)
            ylabel(ax(5),'V^2/Hz','FontSize',20)
            
            
            pause(.0000001)
            delete(l0);
            delete(l1);
            delete(l2);
            delete(l3);
            delete(l4);
            delete(l5);
            delete(l6);
        end
    end
catch
    fclose(ser);
    disp('ca a merde')
    
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
    window = ones(Lax1,1)*hanning(Lscan).';
    wc2=1/mean(window(1,:).^2);            % window correction factor
    datap  = window.*(data- mean(data,2)* ones(1,Lscan));
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
            P12(j,:,:,:)=conj(tempo).*P1(ind_nbsensor~=j,:,:)./Lscan^2/dk*wc2;
            tempo=shiftdim(repmat(squeeze(P11(j,:,:)),[1,1,nb_sensor-1]),2);
            Co12(j,:,:,:)=squeeze(P12(j,:,:,:)).^2./(tempo.*P11(ind_nbsensor~=j,:,:));
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





