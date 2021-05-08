function Accel_freq_mean_deployment(z,MS,channel,ax,titlestr,limitcolor,XTickLabel,ylabelstr,xlabelstr)

    Z=length(z);

    interp_profile=@(x,y) (interp1(x.pr,squeeze(log10(x.Pf(y,:,:))),z));
    makemean=@(x,y,z) (nanmean(reshape(x,[Z,y,z]),3));
    
    ind0=cellfun(@isempty,MS);
    MS=MS(~ind0);

    MSmean=cellfun(@(x) interp_profile(x,channel),MS,'un',0);
    MSmean=cell2mat(MSmean);
    
    F=length(MS{1}.f);
    f=MS{1}.f;
    P=length(MS);
    D=makemean(MSmean,F,P);
    maxZ=median(cellfun(@(x) max(x.pr),MS));
    
    colormap(ax,cmocean('balance'))
    pcolor(ax,f,z,D);shading(ax,'flat')
    caxis(ax,limitcolor)
    grid(ax,'on')
    ylim(ax,[0 maxZ+5])
    xlim(ax,[1/3 180])
    title(ax,titlestr)
    ax.FontSize=20;
    ylabel(ax,ylabelstr,'fontsize',20)
    xlabel(ax,xlabelstr,'fontsize',20)
    set(ax,'XScale','log','YDir','reverse')
    ax.XTick=[1 10 100];
    ax.YTick=floor(linspace(0,maxZ,10));
    if isempty(XTickLabel);ax.XTickLabel=XTickLabel;end;
end
