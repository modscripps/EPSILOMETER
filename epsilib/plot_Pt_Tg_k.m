function [] = plot_Pt_Tg_k(scan)

% Get colors
cols = mod_epsilometer_define_colors;

% Smooth Tdiff wavenumber spectra
smTG1 = smoothdata(scan.Pt_Tg_k.t1,'movmean',10);
smTG2 = smoothdata(scan.Pt_Tg_k.t2,'movmean',10);

% Get Batchelor spectrum for each combination of epsilon/chi
[kbatch11,Pbatch11] = batchelor(scan.epsilon.s1,scan.chi.t1, ...
    scan.kvis,scan.ktemp);
[kbatch21,Pbatch21] = batchelor(scan.epsilon.s2,scan.chi.t1, ...
    scan.kvis,scan.ktemp);

[kbatch12,Pbatch12] = batchelor(scan.epsilon.s1,scan.chi.t2, ...
    scan.kvis,scan.ktemp);
[kbatch22,Pbatch22] = batchelor(scan.epsilon.s2,scan.chi.t2, ...
    scan.kvis,scan.ktemp);

% Get FPO7 noise
Meta_Data = scan.Meta_Data;
dTdV=[Meta_Data.AFE.t1.cal Meta_Data.AFE.t2.cal];

logf=log10(scan.f);
FPO7noise = scan.FPO7noise;
n0=FPO7noise.n0; n1=FPO7noise.n1; n2=FPO7noise.n2; n3=FPO7noise.n3;
tnoise=10.^(n0+n1.*logf+n2.*logf.^2+n3.*logf.^3);

h_freq=get_filters_SOM(Meta_Data,scan.f);
k_noise=scan.f./scan.w;
noise_t=tnoise.*dTdV(1).^2./h_freq.FPO7(scan.w);
tnoise_k= (2*pi*k_noise).^2 .* noise_t.*scan.w;        % T1_k spec  as function of k


%% Make the plot
p9(1) = loglog(scan.k,scan.Pt_Tg_k.t1,':','color',cols.t1,'linewidth',2);
hold on
p9(2) = loglog(scan.k,smTG1,'color',cols.t1,'linewidth',3);
p9(3) = loglog(scan.k,scan.Pt_Tg_k.t2,':','color',cols.t2);
p9(4) = loglog(scan.k,smTG2,'color',cols.t2);

% Add Batchelor spectra
p9(5) = loglog(kbatch11,Pbatch11,'Color',cols.batch11);
p9(6) = loglog(kbatch12,Pbatch12,'Color',cols.batch12);
p9(7) = loglog(kbatch21,Pbatch21,'Color',cols.batch21);
p9(8) = loglog(kbatch22,Pbatch22,'Color',cols.batch22);

% Add noise
p9(9) = loglog(k_noise,tnoise_k,'k:');

% Add kc
indkc=find(scan.k>scan.kc.t1,1,'first');
p9(10) = scatter(scan.k(indkc),smTG1(indkc),'filled','d','sizedata',300,'MarkerEdgeColor','k','markerfacecolor',cols.t1,'linewidth',2);
indkc=find(scan.k>scan.kc.t2,1,'first');
p9(11) = scatter(scan.k(indkc),smTG2(indkc),'filled','p','sizedata',450,'MarkerEdgeColor','k','markerfacecolor',cols.t2,'linewidth',2);

legend('t1','t1smooth','t2','t2smooth',...
    'batch11','batch12','batch21','batch22',...
    'T-noise','t1_{cutoff}','t2_{cutoff}',...
    'location','southwest','numcolumns',3);
xlim([6e-1 400])
ylim([1e-10 1e-1])
grid on
xlabel('k (cpm)')
ylabel('\phi^2_{TG} (C^2 m^{-2} / cpm)')