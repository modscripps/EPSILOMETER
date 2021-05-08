% check TF quality flag


prid=2
Hnoise=interp1(fH,H,MS1{prid}.f);
ind_pump=147
close all
for n=1:MS1{prid}.nbscan
    TFa3=Hnoise.*squeeze(MS1{prid}.Pf(end,n,:))';
    A=smoothdata(squeeze(MS1{prid}.Pf(3,n,:)),'movmean',5);
    B=smoothdata(2*TFa3,'movmean',5);
    B1=B*A(147)/B(147);
    loglog(MS1{10}.f,A,'b');
    hold on
%     loglog(MS1{10}.f,smoothdata(squeeze(MS1{10}.Pf(end,n,:))),'r');
    loglog(MS1{10}.f,B,'k');
    loglog(MS1{10}.f,B1,'m');
    legend('s1','Hnoise')
    title(sprintf('%i over %i',n,MS1{prid}.nbscan))
    ylim([1e-13 1e-2])
    grid on
    qc_flag=log10(A./B.');
    qc_flag=nanmean(qc_flag(MS1{prid}.f>2 & MS1{prid}.f<35));
    qc_flag1=log10(A./B1.');
    qc_flag1=nanmean(qc_flag1(MS1{prid}.f>2 & MS1{prid}.f<35));
    plot([2 2],[1e-13 1e-2],'k--')
    plot([35 35],[1e-13 1e-2],'k--')
    text(.2,1e-4,sprintf('%1.2f',qc_flag),'fontsize',15)
    text(.6,1e-4,sprintf('%1.2f',qc_flag1),'backgroundcolor','m','fontsize',15)
    pause
    clf
end