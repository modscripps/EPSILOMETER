function MS=fitpanchev(Meta_Data,MS,dsp)

%%load('/Volumes/DataDrive/ArrowII_epsi/WW/d1/L1/Turbulence_Profiles.mat');
if nargin<3
    dsp=0;
end
if dsp==1
    v = VideoWriter(sprintf('%s_fitpanchev%s',Meta_Data.deployment,'.avi'));
    v.FrameRate=5;
    open(v)
end

epsilonrange=logspace(-11,-5,100);
for id_profile=1:length(MS)
    fprintf('%i over %i\n',id_profile,length(MS));
    for id_scan=1:MS{id_profile}.nbscan
        clear kpan;clear Ppan;
        for e=1:length(epsilonrange)
            [kpan{e},Ppan{e}] = panchev(epsilonrange(e),MS{id_profile}.kvis(id_scan));
        end
        kc=MS{id_profile}.kc(id_scan,1);
        kall=kpan{1}(1):nanmean(diff(kpan{1})):kc;
        Ppan_array=cell2mat(cellfun(@(x,y) interp1(x,y,kall),kpan,Ppan,'un',0).');
        
        obs1=interp1(MS{id_profile}.k,smoothdata(squeeze(MS{id_profile}.Pshearco_k(id_scan,:,1)),'movmean',10),kall);
        obs2=interp1(MS{id_profile}.k,smoothdata(squeeze(MS{id_profile}.Pshearco_k(id_scan,:,2)),'movmean',10),kall);
        fitPan1=abs(log10(Ppan_array./obs1));
        fitPan2=abs(log10(Ppan_array./obs2));
        for i=1:length(epsilonrange)
            indfit1(i)=length(find(fitPan1(i,:)>-5e-3 & fitPan1(i,:)<5e-3));
            indfit2(i)=length(find(fitPan2(i,:)>-5e-3 & fitPan2(i,:)<5e-3));
            %ind1{i}=find(fitPan1(i,:)>-5e-3 & fitPan1(i,:)<5e-3);
        end
        smindfit1=smoothdata(indfit1,'movmean',10);
        smindfit2=smoothdata(indfit2,'movmean',10);
        [pks1,locs1]= findpeaks(smindfit1);
        [pks2,locs2]= findpeaks(smindfit2);
        
        if ~isempty(locs1)
            new_epsilon{id_profile}(id_scan,1)= epsilonrange(locs1(1));
            
            if dsp==1
%                 figure(1)
%                 plot(log10(epsilonrange),indfit1)
%                 hold on
%                 plot(log10(epsilonrange),smindfit,'m')
%                 hold off
                
                loglog(MS{id_profile}.k,smoothdata(squeeze(MS{id_profile}.Pshearco_k(id_scan,:,1)),10),'b');
                hold on
                %indlocs1=find(cellfun(@isempty,ind1(locs1))==0,1,'first');
                %loglog(kall(min([ind1{locs1(indlocs1)}]):max([ind1{locs1(indlocs1)}])),obs1(min([ind1{locs1(indlocs1)}]):max([ind1{locs1(indlocs1)}])),'g');
                loglog(kpan{locs1(1)},Ppan{locs1(1)},'m')
                loglog(MS{id_profile}.k,squeeze(MS{id_profile}.Ppan(id_scan,:,1)),'k')
                title(sprintf('scan%i-New %1.3e eps/Old %1.3e eps',id_scan,epsilonrange(locs1(1)),MS{id_profile}.epsilon_old(id_scan,1)))
                legend('data','new panchev','old panchev')
                hold off
                pause(0.1)
                xlim([1e-1 1e3])
                xlabel('cpm','fontsize',20)
                ylabel('shear^2/cpm','fontsize',20)
                % movie stuff
                frame=getframe(gcf);
                writeVideo(v,frame)
                cla
            end
        else
            new_epsilon{id_profile}(id_scan,1)= nan;
        end
        if ~isempty(locs2)
            new_epsilon{id_profile}(id_scan,2)= epsilonrange(locs2(1));
        else
            new_epsilon{id_profile}(id_scan,2)= nan;
        end
    end
end

%%



if dsp==1
    close(v)
else
    for id_profile=1:length(MS)
        fprintf('%i over %i\n',id_profile,length(MS));
        MS{id_profile}.epsilon_fit=new_epsilon{id_profile};
    end
    %save(fullfile(Meta_Data.L1path,'Turbulence_Profiles.mat'),'MS','-v7.3')
end

end

