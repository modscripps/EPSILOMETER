function plot_binned_epsilon_sanity_profile(Epsilon_class,title_string,ax1,ax2)

cmap = parula(length(Epsilon_class.bin)+1);

pp = loglog(ax1,Epsilon_class.kpan.',Epsilon_class.Ppan.','color',[1 1 1]*.7,'linewidth',1);
ll=pp;

mink=Epsilon_class.k(find(~isnan(nansum(Epsilon_class.mPshear1,1)),1,'first'));
maxk=Epsilon_class.k(find(~isnan(nansum(Epsilon_class.mPshear1,1)),1,'last'));
minP=min(Epsilon_class.mPshear1(:));
%maxP=max(Epsilon_class.mPshear1(:));
maxP=max(Epsilon_class.Ppan(:)).*10;

hold(ax1,'on')
ind_bin=find(Epsilon_class.nbin1>0);
%ind_bin=[];%trick for SP1810
ind_bin=ind_bin(1:end-2);
for i=ind_bin
    indk=find(~isnan(Epsilon_class.Ppan(i,:)),1,'first');
    if Epsilon_class.kpan(i,indk)<=mink
        text(mink, ...
            Epsilon_class.Ppan(i,indk),...
            sprintf('10^{%1.1f}',log10(Epsilon_class.bin(i))),...
            'fontsize',15,'Parent',ax1)
    else
        text(Epsilon_class.kpan(i,indk), ...
            Epsilon_class.Ppan(i,indk),...
            sprintf('10^{%1.1f}',log10(Epsilon_class.bin(i))),...
            'fontsize',15,'Parent',ax1)
    end
    ll(i)=loglog(ax1,Epsilon_class.k,Epsilon_class.mPshear1(i,:),'Color',cmap(i,:),'linewidth',2);
end
text(10,5.*max(Epsilon_class.Ppan(:)),'Shear1','fontsize',15,'Parent',ax1)


xlim(ax1,[mink maxk])
ylim(ax1,[minP maxP])
grid(ax1,'on')
%xlabel(ax1,'k [cpm]')
ylabel(ax1,'[s^{-2} / cpm]')
set(ax1,'fontsize',15)
set(ax1,'Xticklabel','')

%title(ax1,[title_string ' Shear1'])

nanind=find(nansum(Epsilon_class.mPshear1,2)>0);
%legend_string=arrayfun(@(x) sprintf('log_{10}(\\epsilon)~%2.1f',x),log10(Epsilon_class.bin),'un',0);
%hl=gridLegend(ll(nanind),2,legend_string{nanind},'location','northwest');
%set(hl,'fontsize',10)

cmap = parula(length(Epsilon_class.bin)+1);

pp = loglog(ax2,Epsilon_class.kpan.',Epsilon_class.Ppan.','color',[1 1 1]*.7,'linewidth',1);
ll=pp;

mink=Epsilon_class.k(find(~isnan(nansum(Epsilon_class.mPshear2,1)),1,'first'));
maxk=Epsilon_class.k(find(~isnan(nansum(Epsilon_class.mPshear2,1)),1,'last'));
minP=min(Epsilon_class.mPshear2(:));
%maxP=max(Epsilon_class.mPshear2(:));
maxP=max(Epsilon_class.Ppan(:)).*10;

hold(ax2,'on')
ind_bin=find(Epsilon_class.nbin2>0);
ind_bin=ind_bin(1:end-2);
for i=ind_bin
    indk=find(~isnan(Epsilon_class.Ppan(i,:)),1,'first');
    if Epsilon_class.kpan(i,indk)<=mink
        text(mink, ...
            Epsilon_class.Ppan(i,indk),...
            sprintf('10^{%1.1f}',log10(Epsilon_class.bin(i))),...
            'fontsize',15,'Parent',ax2)
    else
        text(Epsilon_class.kpan(i,indk), ...
            Epsilon_class.Ppan(i,indk),...
            sprintf('10^{%1.1f}',log10(Epsilon_class.bin(i))),...
            'fontsize',15,'Parent',ax2)
    end
    ll(i)=loglog(ax2,Epsilon_class.k,Epsilon_class.mPshear2(i,:),'Color',cmap(i,:),'linewidth',2);
end
text(10,5.*max(Epsilon_class.Ppan(:)),'Shear2','fontsize',15,'Parent',ax2)

xlim(ax2,[mink maxk])
ylim(ax2,[minP maxP])
grid(ax2,'on')
set(ax2,'Xscale','log','Yscale','log')
xlabel(ax2,'k [cpm]')
ylabel(ax2,'[s^{-2} / cpm]')
set(ax2,'fontsize',15)
%title(ax2,[title_string ' Shear2'])

% nanind=find(nansum(Epsilon_class.mPshear2,2)>0);
% legend_string=arrayfun(@(x) sprintf('log_{10}(\\epsilon)~%2.1f',x),log10(Epsilon_class.bin),'un',0);
% hl=gridLegend(ll(nanind),2,legend_string{nanind},'location','northwest');
% set(hl,'fontsize',10)

