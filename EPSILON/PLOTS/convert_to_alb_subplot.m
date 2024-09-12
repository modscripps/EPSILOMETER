function ax=convert_to_alb_subplot(ax,Labelfontsize,axefontsize)

for a=1:length(ax)
    ax(a).FontSize=axefontsize;
    ax(a).XLabel.FontSize=Labelfontsize;
    ax(a).FontName='Times New Roman';
    ax(a).XLabel.FontName='Times New Roman';
end

