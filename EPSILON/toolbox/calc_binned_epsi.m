function Epsilon_class=calc_binned_epsi(MS,epsi_bin)


%  compute averaged shear spectra classed by epsilon values
%
%  input: 
% . MS: structure with the chi values and the temperature gradiant spectra
% . epsi_bin: bin of epsilon values
%
%  Created by Arnaud Le Boyer on 7/28/18.

    %% create epsi_bin or include min and max epsi value to epsi_bin
    if nargin<2
        epsi_bin=10.^(-10:.5:-4.5);
    end
    
    % trick for conditional anonymous function
    iif = @(varargin) varargin{2 * find([varargin{1:2:end}], 1, 'first')}();

    
    % change struct to cell TODO check if it would be worth it to do it for
    %  the whole process
    if iscell(MS)
        S_MS=[MS{:}];
    else
        S_MS=MS;
    end
        
    %% create common k(cpm) axis
    maxk = min([max([S_MS.k]),300]);
    mink = min([S_MS.k]);
    if mink==0
        mink=sort([S_MS.k]);
        mink=mink(find(mink>0,1,'first'));
    end
    dk   = max(mink,.1);
    
    k    = mink:dk:maxk;
    
    % project Epsilon fields onto the common k and omega axis 
    Psheark=cellfun(@(x) shiftdim(interp1(x.k,shiftdim(x.Pshear_k,1),k),1) ,MS,'un',0);
    Psheark_co=cellfun(@(x) shiftdim(interp1(x.k,shiftdim(x.Pshearco_k,1),k),1) ,MS,'un',0);
%     Psheark_eof=cellfun(@(x) shiftdim(interp1(x.k,shiftdim(x.Psheareof_k,1),k),1) ,MS,'un',0);
    Psheark=[Psheark{:}];
    Psheark_co=[Psheark_co{:}];
%     Psheark_eof=[Psheark_eof{:}];
    
    epsilon=cat(1,S_MS.epsilon);
    epsilon_co=cat(1,S_MS.epsilon_co);
%     epsilon_eof=cat(1,S_MS.epsilon_eof);
    w=cat(2,S_MS.w);
    kvis=real(cat(1,S_MS.kvis));
    
%     index1=arrayfun(@(x) find(epsilon(:,1)>=x-.5*x & epsilon(:,1)<x+.5*x & w.'>.3 & w.'<.7),epsi_bin,'un',0);
%     index2=arrayfun(@(x) find(epsilon(:,2)>=x-.5*x & epsilon(:,2)<x+.5*x & w.'>.3 & w.'<.7),epsi_bin,'un',0);
%     index1_co=arrayfun(@(x) find(epsilon_co(:,1)>=x-.5*x & epsilon_co(:,1)<x+.5*x & w.'>.3 & w.'<.7),epsi_bin,'un',0);
%     index2_co=arrayfun(@(x) find(epsilon_co(:,2)>=x-.5*x & epsilon_co(:,2)<x+.5*x & w.'>.3 & w.'<.7),epsi_bin,'un',0);
    index1=arrayfun(@(x) find(epsilon(:,1)>=x-.5*x & epsilon(:,1)<x+.5*x),epsi_bin,'un',0);
    index2=arrayfun(@(x) find(epsilon(:,2)>=x-.5*x & epsilon(:,2)<x+.5*x),epsi_bin,'un',0);
    index1_co=arrayfun(@(x) find(epsilon_co(:,1)>=x-.5*x & epsilon_co(:,1)<x+.5*x),epsi_bin,'un',0);
    index2_co=arrayfun(@(x) find(epsilon_co(:,2)>=x-.5*x & epsilon_co(:,2)<x+.5*x),epsi_bin,'un',0);
%     index1_eof=arrayfun(@(x) find(epsilon_eof(:,1)>=x-.5*x & epsilon_eof(:,1)<x+.5*x),epsi_bin,'un',0);
%     index2_eof=arrayfun(@(x) find(epsilon_eof(:,2)>=x-.5*x & epsilon_eof(:,2)<x+.5*x),epsi_bin,'un',0);
    %% if thehre is only 1 epsilon value in the bin, empty the bin.
    index1=cellfun(@(x) iif(length(x)<=1,[],length(x)>1,x),index1,'un',0);
    index2=cellfun(@(x) iif(length(x)<=1,[],length(x)>1,x),index2,'un',0);
    index1_co=cellfun(@(x) iif(length(x)<=1,[],length(x)>1,x),index1_co,'un',0);
    index2_co=cellfun(@(x) iif(length(x)<=1,[],length(x)>1,x),index2_co,'un',0);
%     index1_eof=cellfun(@(x) iif(length(x)<=1,[],length(x)>1,x),index1_eof,'un',0);
%     index2_eof=cellfun(@(x) iif(length(x)<=1,[],length(x)>1,x),index2_eof,'un',0);
    
    %% define Epsilon class
    % epislon values
    % get and average the spectra per epsilon values, 
    % get the number of spectra per epsilon values,
    % get the median kinematic viscosity per epsilon values
    % compute the panchev spectra for epsilon values
    Epsilon_class=struct();
    Epsilon_class.k=k;
    Epsilon_class.bin=epsi_bin;
    
    Epsilon_class.nbin1=cellfun(@length,index1);
    Epsilon_class.nbin2=cellfun(@length,index2);
    Epsilon_class.nbin1_co=cellfun(@length,index1_co);
    Epsilon_class.nbin2_co=cellfun(@length,index2_co);
%     Epsilon_class.nbin1_eof=cellfun(@length,index1_eof);
%     Epsilon_class.nbin2_eof=cellfun(@length,index2_eof);
    
    Epsilon_class.Pshear1=cellfun(@(x) squeeze(Psheark(1,x,:)),index1,'un',0);
    Epsilon_class.Pshear2=cellfun(@(x) squeeze(Psheark(2,x,:)),index2,'un',0);
    Epsilon_class.Pshear1_co=cellfun(@(x) squeeze(Psheark_co(1,x,:)),index1_co,'un',0);
    Epsilon_class.Pshear2_co=cellfun(@(x) squeeze(Psheark_co(2,x,:)),index2_co,'un',0);
%     Epsilon_class.Pshear1_eof=cellfun(@(x) squeeze(Psheark_eof(1,x,:)),index1_eof,'un',0);
%     Epsilon_class.Pshear2_eof=cellfun(@(x) squeeze(Psheark_eof(2,x,:)),index2_eof,'un',0);
    
    Epsilon_class.mPshear1=cell2mat(cellfun(@(x) nanmean(x,1),Epsilon_class.Pshear1,'un',0).');
    Epsilon_class.mPshear2=cell2mat(cellfun(@(x) nanmean(x,1),Epsilon_class.Pshear2,'un',0).');
    Epsilon_class.mPshear1_co=cell2mat(cellfun(@(x) nanmean(x,1),Epsilon_class.Pshear1_co,'un',0).');
    Epsilon_class.mPshear2_co=cell2mat(cellfun(@(x) nanmean(x,1),Epsilon_class.Pshear2_co,'un',0).');
%     Epsilon_class.mPshear1_eof=cell2mat(cellfun(@(x) nanmean(x,1),Epsilon_class.Pshear1_eof,'un',0).');
%     Epsilon_class.mPshear2_eof=cell2mat(cellfun(@(x) nanmean(x,1),Epsilon_class.Pshear2_eof,'un',0).');
    
    Epsilon_class.kvis=cellfun(@(x) median(kvis(x).'),index1,'un',0);
    Epsilon_class.kvis=cellfun(@(x) iif(isnan(x),nanmedian([Epsilon_class.kvis{:}].'),~isnan(x),x), ...
                   Epsilon_class.kvis,'un',0);
   
    [kpan,Ppan] = cellfun(@(x,y) panchev(x,y),num2cell(epsi_bin),Epsilon_class.kvis,'un',0);
    Epsilon_class.kpan = cell2mat(kpan).';
    Epsilon_class.Ppan = cell2mat(Ppan).';
    
end

