function Chi_class=calc_binned_chi(MS,epsi_bin)

%  compute averaged temperature gradiant spectra calssed by chi values
%
%  input: 
% . MS: structure with the chi values and the temperature gradiant spectra
% . epsi_bin: bin a chi and epsilon values
%
%  Created by Arnaud Le Boyer on 7/28/18.



    %% create epsi_bin or include min and max epsi value to epsi_bin
    if nargin<2
        chi_bin =10.^(-10:.5:-4.5);
        epsi_bin=10.^(-10:.5:-4.5);
    end
    

    
    % change struct to cell TODO check if it would be worth it to do it for
    %  the whole process
    if iscell(MS)
        S_MS=[MS{:}];
    end
        
    %% create common k(cpm) axis
    ktot=sort(unique([S_MS.k]));
    maxk = min([max([S_MS.k]),300]);
    mink = min([S_MS.k]);
    dk   = max(ktot(2),.1);
    k    = mink:dk:maxk;
    
    % project Epsilon fields onto the common k and omega axis 
    PphiT_k=cellfun(@(x) shiftdim(interp1(x.k,shiftdim(x.PphiT_k,1),k),1) ,MS,'un',0);
    PphiT_k=[PphiT_k{:}];
    
    chi=cat(1,S_MS.chi);
    epsilon=cat(1,S_MS.epsilon_co);
    kvis=real(cat(1,S_MS.kvis));
    ktemp=real(cat(1,S_MS.ktemp));

    index11=cell(length(chi_bin),length(epsi_bin));
    index12=cell(length(chi_bin),length(epsi_bin));
    index21=cell(length(chi_bin),length(epsi_bin));
    index22=cell(length(chi_bin),length(epsi_bin));
    for i=1:length(chi_bin)
        x=chi_bin(i);
        for j=1:length(epsi_bin)
            y=epsi_bin(j);
            ind=find(chi(:,1)>=x-.5*x & chi(:,1)<x+.5*x & ...
                         epsilon(:,1)>=y-.5*y & epsilon(:,1)<y+.5*y);
            index11{i,j}=ind;
            ind=find(chi(:,1)>=x-.5*x & chi(:,1)<x+.5*x & ...
                         epsilon(:,2)>=y-.5*y & epsilon(:,2)<y+.5*y);
            index12{i,j}=ind;
            ind=find(chi(:,2)>=x-.5*x & chi(:,2)<x+.5*x & ...
                         epsilon(:,1)>=y-.5*y & epsilon(:,1)<y+.5*y);
            index21{i,j}=ind;
            ind=find(chi(:,2)>=x-.5*x & chi(:,2)<x+.5*x & ...
                         epsilon(:,2)>=y-.5*y & epsilon(:,2)<y+.5*y);
            index22{i,j}=ind;
        end
    end
    %% if thehre is only 1 epsilon value in the bin, empty the bin.
    % trick for conditional anonymous function
    iif = @(varargin) varargin{2 * find([varargin{1:2:end}], 1, 'first')}();
    
    index11=cellfun(@(x) iif(length(x)<=1,[],length(x)>1,x),index11,'un',0);
    index12=cellfun(@(x) iif(length(x)<=1,[],length(x)>1,x),index12,'un',0);
    index21=cellfun(@(x) iif(length(x)<=1,[],length(x)>1,x),index21,'un',0);
    index22=cellfun(@(x) iif(length(x)<=1,[],length(x)>1,x),index22,'un',0);
    
    Chi_class=struct();
    Chi_class.k=k;
    Chi_class.epsibin=epsi_bin;
    Chi_class.chibin=chi_bin;
    Chi_class.nbin11=cellfun(@length,index11);
    Chi_class.nbin12=cellfun(@length,index12);
    Chi_class.nbin21=cellfun(@length,index21);
    Chi_class.nbin22=cellfun(@length,index22);

    Chi_class.PphiT11=cellfun(@(x) squeeze(PphiT_k(1,x,:)),index11,'un',0);
    Chi_class.PphiT21=cellfun(@(x) squeeze(PphiT_k(2,x,:)),index21,'un',0);
    Chi_class.PphiT12=cellfun(@(x) squeeze(PphiT_k(1,x,:)),index12,'un',0);
    Chi_class.PphiT22=cellfun(@(x) squeeze(PphiT_k(2,x,:)),index22,'un',0);

    Chi_class.mPphiT11=cellfun(@(x) nanmean(x,1),Chi_class.PphiT11,'un',0);
    Chi_class.mPphiT21=cellfun(@(x) nanmean(x,1),Chi_class.PphiT21,'un',0);
    Chi_class.mPphiT12=cellfun(@(x) nanmean(x,1),Chi_class.PphiT12,'un',0);
    Chi_class.mPphiT22=cellfun(@(x) nanmean(x,1),Chi_class.PphiT22,'un',0);

    Chi_class.mPphiT11=reshape(vertcat(Chi_class.mPphiT11{:}),[length(chi_bin),length(epsi_bin),length(k)]);
    Chi_class.mPphiT12=reshape(vertcat(Chi_class.mPphiT12{:}),[length(chi_bin),length(epsi_bin),length(k)]);
    Chi_class.mPphiT21=reshape(vertcat(Chi_class.mPphiT21{:}),[length(chi_bin),length(epsi_bin),length(k)]);
    Chi_class.mPphiT22=reshape(vertcat(Chi_class.mPphiT22{:}),[length(chi_bin),length(epsi_bin),length(k)]);
    
    Chi_class.kvis11=cellfun(@(x) nanmean(kvis(x)),index11,'un',0);
    Chi_class.kvis12=cellfun(@(x) nanmean(kvis(x)),index12,'un',0);
    Chi_class.kvis21=cellfun(@(x) nanmean(kvis(x)),index21,'un',0);
    Chi_class.kvis22=cellfun(@(x) nanmean(kvis(x)),index22,'un',0);
    
    Chi_class.ktemp11=cellfun(@(x) nanmean(ktemp(x)),index11,'un',0);
    Chi_class.ktemp12=cellfun(@(x) nanmean(ktemp(x)),index12,'un',0);
    Chi_class.ktemp21=cellfun(@(x) nanmean(ktemp(x)),index21,'un',0);
    Chi_class.ktemp22=cellfun(@(x) nanmean(ktemp(x)),index22,'un',0);
    
    
    Chi_class.epsilon11=cellfun(@(x) nanmean(epsilon(x)),index11,'un',0);
    Chi_class.epsilon12=cellfun(@(x) nanmean(epsilon(x)),index12,'un',0);
    Chi_class.epsilon21=cellfun(@(x) nanmean(epsilon(x)),index21,'un',0);
    Chi_class.epsilon22=cellfun(@(x) nanmean(epsilon(x)),index22,'un',0);
    
    Chi_class.kvis11=cellfun(@(x) iif(isnan(x),nanmean([Chi_class.kvis11{:}]),~isnan(x),x), ...
                   Chi_class.kvis11,'un',0);
    Chi_class.kvis12=cellfun(@(x) iif(isnan(x),nanmean([Chi_class.kvis12{:}]),~isnan(x),x), ...
                   Chi_class.kvis12,'un',0);
    Chi_class.kvis21=cellfun(@(x) iif(isnan(x),nanmean([Chi_class.kvis21{:}]),~isnan(x),x), ...
                   Chi_class.kvis21,'un',0);
    Chi_class.kvis22=cellfun(@(x) iif(isnan(x),nanmean([Chi_class.kvis22{:}]),~isnan(x),x), ...
                   Chi_class.kvis22,'un',0);

    Chi_class.ktemp11=cellfun(@(x) iif(isnan(x),nanmean([Chi_class.ktemp11{:}]),~isnan(x),x), ...
                   Chi_class.ktemp11,'un',0);
    Chi_class.ktemp12=cellfun(@(x) iif(isnan(x),nanmean([Chi_class.ktemp12{:}]),~isnan(x),x), ...
                   Chi_class.ktemp12,'un',0);
    Chi_class.ktemp21=cellfun(@(x) iif(isnan(x),nanmean([Chi_class.ktemp21{:}]),~isnan(x),x), ...
                   Chi_class.ktemp21,'un',0);
    Chi_class.ktemp22=cellfun(@(x) iif(isnan(x),nanmean([Chi_class.ktemp22{:}]),~isnan(x),x), ...
                   Chi_class.ktemp22,'un',0);

               
    Chi_class.kbatch=k;

    for i=1:length(chi_bin)
        x=chi_bin(i);
        for j=1:length(epsi_bin)
            y=epsi_bin(j);
            [kbatch,Pbatch] = batchelor(y,x, ...
                               Chi_class.kvis11{i,j},Chi_class.ktemp11{i,j});
            if ~isnan(kbatch)
                Chi_class.Pbatch11(i,j,:)=interp1(kbatch,Pbatch,k);
            else
                Chi_class.Pbatch11(i,j,:)=k.*nan;
            end
            [kbatch,Pbatch] = batchelor(y,x, ...
                               Chi_class.kvis12{i,j},Chi_class.ktemp12{i,j});
            if ~isnan(kbatch)
                Chi_class.Pbatch12(i,j,:)=interp1(kbatch,Pbatch,k);
            else
                Chi_class.Pbatch12(i,j,:)=k.*nan;
            end
            [kbatch,Pbatch] = batchelor(y,x, ...
                               Chi_class.kvis21{i,j},Chi_class.ktemp21{i,j});
            if ~isnan(kbatch)
                Chi_class.Pbatch21(i,j,:)=interp1(kbatch,Pbatch,k);
            else
                Chi_class.Pbatch21(i,j,:)=k.*nan;
            end
            [kbatch,Pbatch] = batchelor(y,x, ...
                               Chi_class.kvis22{i,j},Chi_class.ktemp22{i,j});
            if ~isnan(kbatch)
                Chi_class.Pbatch22(i,j,:)=interp1(kbatch,Pbatch,k);
            else
                Chi_class.Pbatch22(i,j,:)=k.*nan;
            end
        end
    end
            
    

