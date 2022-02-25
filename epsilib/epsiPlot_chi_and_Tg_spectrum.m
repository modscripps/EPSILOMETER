function epsiPlot_chi_and_Tg_spectrum(Profile,depth)

%% Setup axes
figure('units','inches','position',[0 0 17.8 15]);
ax(1) = subtightplot(4,4,[4,8]);
ax(2) = subtightplot(4,4,[1,5]);
ax(3) = subtightplot(4,4,[2:3,6:7]);

ax(4) = subtightplot(4,4,[12,16]);
ax(5) = subtightplot(4,4,[9,13]);
ax(6) = subtightplot(4,4,[10:11,14:15]);

pp = epsiSetup_set_plot_properties;
cols = pp.Colors;
cols.w = [0 0 0];
cols.profInfo = [174 197 227]./255;
cols.scanInfo = [0.9 0.9 0.9];
cols.batch_s1t1 = [0.4902    0.1647    0.4706];
cols.batch_s1t2 = [0.6784    0.1529    0.6431];
cols.batch_s2t1 = [0.8353    0.1059    0.7922];
cols.batch_s2t2 = [0.9679    0.4079    0.6317];
cols.panchev1 = [0.4902    0.1647    0.4706];
cols.panchev2 = [0.8353    0.1059    0.7922];

%% Scan data
[~,k] = min(abs(Profile.pr - depth));
scan = get_scan_spectra(Profile,k);

if isfield(scan,'Pt_Tg_k')
% Smooth Tdiff wavenumber spectra
smTG1 = smoothdata(scan.Pt_Tg_k.t1,'movmean',10);
smTG2 = smoothdata(scan.Pt_Tg_k.t2,'movmean',10);

%% Profile data
Meta_Data = Profile.Meta_Data;
dTdV=[Meta_Data.AFE.t1.cal Meta_Data.AFE.t2.cal];

Profile.epsi.P = interp1(Profile.ctd.time_s,Profile.ctd.P,Profile.epsi.time_s);
prInScan = Profile.epsi.P(scan.ind_scan);
prRange = [nanmin(prInScan) nanmax(prInScan)];

% noise floor
logf=log10(Profile.f);
h_freq=get_filters_SOM(Meta_Data,Profile.f);
switch Meta_Data.AFE.temp_circuit
    case 'Tdiff'
        FPO7noise=load(fullfile(Meta_Data.paths.calibration,'FPO7_noise.mat'),'n0','n1','n2','n3');
    otherwise
        FPO7noise=load(fullfile(Meta_Data.paths.calibration,'FPO7_notdiffnoise.mat'),'n0','n1','n2','n3');
end
n0=FPO7noise.n0; n1=FPO7noise.n1; n2=FPO7noise.n2; n3=FPO7noise.n3;
tnoise=10.^(n0+n1.*logf+n2.*logf.^2+n3.*logf.^3);

%% Get FPO7 noise
k_noise=Profile.f./Profile.w(k);
noise_t=tnoise.*dTdV(1).^2./h_freq.FPO7(Profile.w(k));
tnoise_k= (2*pi*k_noise).^2 .* noise_t.*Profile.w(k);        % T1_k spec  as function of k

% Get Batchelor spectrum for each combination of epsilon/chi
[kbatch_s1t1,Pbatch_s1t1] = batchelor(scan.epsilon.s1,scan.chi.t1, ...
    scan.kvis,scan.ktemp);
[kbatch_s2t1,Pbatch_s2t1] = batchelor(scan.epsilon.s2,scan.chi.t1, ...
    scan.kvis,scan.ktemp);

[kbatch_s1t2,Pbatch_s1t2] = batchelor(scan.epsilon.s1,scan.chi.t2, ...
    scan.kvis,scan.ktemp);
[kbatch_s2t2,Pbatch_s2t2] = batchelor(scan.epsilon.s2,scan.chi.t2, ...
    scan.kvis,scan.ktemp);

%% Plot chi profiles
axes(ax(1))
plot(Profile.chi(:,1),Profile.pr,'color',cols.t1,'displayname','\chi_1');
ax(1).XScale = 'log';
ax(1).YLabel.String = 'Depth (m)';
legend('autoupdate','off','location','northwest')
yline(scan.pr,'linewidth',2,'linestyle','--')
ax(1).XTick = logspace(-11,-2,10);

axes(ax(4))
plot(Profile.chi(:,2),Profile.pr,'color',cols.t2,'displayname','\chi_2');
ax(4).XScale = 'log';
ax(4).YLabel.String = 'Depth (m)';
legend('autoupdate','off','location','northwest')
yline(scan.pr,'linewidth',2,'linestyle','--')
ax(4).XTick = logspace(-11,-2,10);

%% Plot t1_volt and t2_volt
indScan = Profile.ind_range_epsi(k,1):Profile.ind_range_epsi(k,2);
axes(ax(2))
plot(scan.t1_volt,Profile.epsi.dnum(scan.ind_scan),'color',cols.t1,'displayname','t1\_volt');
datetick(ax(2),'y','MM:SS')
legend('autoupdate','off','location','northwest')
yline(scan.dnum,'linewidth',2,'linestyle','--')

axes(ax(5))
plot(scan.t2_volt,Profile.epsi.dnum(scan.ind_scan),'color',cols.t2,'displayname','t2\_volt');
datetick(ax(5),'y','MM:SS')
legend('autoupdate','off','location','northwest')
yline(scan.dnum,'linewidth',2,'linestyle','--')

[ax([1:2,4:5]).YDir] = deal('reverse');

%% Plot Tdiff spectra (vs wavenumber)
axes(ax(3))
loglog(scan.k,scan.Pt_Tg_k.t1,':','color',cols.t1,'linewidth',2,'displayname','t1');
hold on
loglog(scan.k,smTG1,'color',cols.t1,'linewidth',3,'displayname','t1\_smooth');
% Add Batchelor spectra
loglog(kbatch_s1t1,Pbatch_s1t1,'Color',cols.batch_s1t1,'displayname','Batchelor s1t1');
loglog(kbatch_s2t1,Pbatch_s2t1,'Color',cols.batch_s2t1,'displayname','Batchelor s2t1');
% Add noise
loglog(k_noise,tnoise_k,'k:','displayname','tnoise');
% Add kc
indkc=find(scan.k>scan.kc.t1,1,'first');
scatter(scan.k(indkc),smTG1(indkc),'filled','p','sizedata',450,...
    'MarkerEdgeColor','k','markerfacecolor','y','linewidth',2,'displayname','t1_{cutoff}');
legend('location','northwest')
grid on

axes(ax(6))
loglog(scan.k,scan.Pt_Tg_k.t2,':','color',cols.t2,'displayname','t2');
hold on
loglog(scan.k,smTG2,'color',cols.t2,'displayname','t2\_smooth');
% Add Batchelor spectra
loglog(kbatch_s1t2,Pbatch_s1t2,'Color',cols.batch_s1t2,'displayname','Batchelor s1t2');
loglog(kbatch_s2t2,Pbatch_s2t2,'Color',cols.batch_s2t2,'displayname','Batchelor s2t2');
% Add noise
loglog(k_noise,tnoise_k,'k:','displayname','tnoise');
% Add kc
indkc=find(scan.k>scan.kc.t2,1,'first');
scatter(scan.k(indkc),smTG2(indkc),'filled','p','sizedata',450,...
    'MarkerEdgeColor','k','markerfacecolor','y','linewidth',2,'displayname','t2_{cutoff}');
legend('location','northwest')
grid on
end