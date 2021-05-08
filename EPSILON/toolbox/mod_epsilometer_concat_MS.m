function [MS,min_epsi,max_epsi]=mod_epsilometer_concat_MS(Meta_Data)


listfile=dir(fullfile(Meta_Data.L1path,'Turbulence_Profiles*.mat'));
listfilename=natsort({listfile.name});
count=1;
channels=Meta_Data.PROCESS.channels;
for f=1:length(listfilename)
    load(fullfile(listfile(f).folder,listfilename{f}),'nb_profile_perfile')
    for p=1:nb_profile_perfile
        load(fullfile(listfile(f).folder,listfilename{f}),sprintf('Profile%03i',count))
        fprintf('%s:%s\r\n',fullfile(listfile(f).folder,listfilename{f}),sprintf('Profile%03i',count));
        eval(sprintf('Profile=Profile%03i;',count));
        Fnames=fieldnames(Profile);
        for n=1:length(Fnames)
            wh_field=Fnames{n};
            switch wh_field
                case 'epsilon'
                    MS{count}.epsilon=Profile.epsilon;
                    minepsi1{count}=nanmin(log10(Profile.epsilon(:,1)));
                    minepsi2{count}=nanmin(log10(Profile.epsilon(:,2)));
                    maxepsi1{count}=nanmax(log10(Profile.epsilon(:,1)));
                    maxepsi2{count}=nanmax(log10(Profile.epsilon(:,2)));
                case 'epsilonTF'
                    MS{count}.epsilonTF=Profile.epsilonTF;
                    minepsi1TF{count}=nanmin(log10(Profile.epsilonTF(:,1)));
                    minepsi2TF{count}=nanmin(log10(Profile.epsilonTF(:,2)));
                    maxepsi1TF{count}=nanmax(log10(Profile.epsilonTF(:,1)));
                    maxepsi2TF{count}=nanmax(log10(Profile.epsilonTF(:,2)));
                case 'Pc1c2'
                    for c=1:length(channels)
                        wh_channel=channels{c};
                        MS{count}.(sprintf('P%s',wh_channel))=Profile.Pc1c2.(wh_channel);
                    end
                case 'Cc1c2'
                    FnamesCc1c2=fieldnames(Profile.Cc1c2);
                    for c=1:length(FnamesCc1c2)
                        wh_channel=FnamesCc1c2{c};
                        MS{count}.(sprintf('C%s',wh_channel))=Profile.Cc1c2.(wh_channel);
                    end
                otherwise
                    MS{count}.(wh_field)=Profile.(wh_field);
                    
            end
        end
        clear Profile;
        eval(sprintf('clear Profile%03i;',count));
        count=count+1;
    end
end
if exist('minepsi1TF','var')
    min_epsi=[nanmean([minepsi1TF{:}]) nanmean([minepsi2TF{:}])];
    max_epsi=[nanmean([maxepsi1TF{:}]) nanmean([maxepsi2TF{:}])];
else
    min_epsi=[nanmean([minepsi1{:}]) nanmean([minepsi2{:}])];
    max_epsi=[nanmean([maxepsi1{:}]) nanmean([maxepsi2{:}])];
end

