function [ax,ax3,p5,p6,p7,p8,p9,p10,f] = plot_profile_scan_spectra2(Meta_Data,Profile_or_profNum,scanNum)
% [p5,p6,p7,p8,p9,p10,f] = plot_profile_scan_spectra(Meta_Data,Profile_or_profNum,scanNum)
%
% p5,p6,p7,p8,p9,p10,f are outputs needed to run this script when making a movie with
% movie_profile_scan_spectra.m

%% Get Profile from Profile_or_profNum
% ----------------------------------------------------------

if isnumeric(Profile_or_profNum) && ~isstruct(Profile_or_profNum)
    profNum = Profile_or_profNum;
    load(fullfile(Meta_Data.L1path,sprintf('Profile%03.0f',profNum)));
    eval(['Profile = ' sprintf('Profile%03.0f',profNum) ';']);
elseif isstruct(Profile_or_profNum)
    Profile = Profile_or_profNum;
    profNum = Profile.profNum;
end

%% Set up figure
% ----------------------------------------------------------

close all
fig = figure;

% Set figure size based on screen size
defaultFigWidth = 1680;
defaultFigHeight = 954;
screenSize = get(0,'screensize');
mult = round(min([screenSize(3)/defaultFigWidth,screenSize(4)/defaultFigHeight]),2);
set(gcf,'Units','pixels','Position',[1 1 defaultFigWidth*mult defaultFigHeight*mult]);

% Get font size based on figure size
axFontSize = floor(16*mult);
infoFontSize = floor(10*mult);

ax(1) = axes('Position',[0.0400    0.5739    0.095    0.4011]); %chi
ax(2) = axes('Position',[0.1555    0.5739    0.095    0.4011]); %epsilon
ax(3) = axes('Position',[0.2700    0.5739    0.095  0.4011]); %T/S
ax(4) = axes('Position',[0.3850    0.5739    0.095    0.4011]); %w
ax(5) = axes('Position',[0.52     0.780     0.215      0.195]); %shear spectra
ax(6) = axes('Position',[0.745    0.780    0.215    0.195]); %coh w/ s1
ax(7) = axes('Position',[0.5200    0.5739    0.215    0.19500]);%acc spectra
ax(8) = axes('Position',[0.745    0.5739    0.215    0.195]); %coh w/ s2
ax(9) = axes('Position',[0.0400    0.14    0.4400    0.3781]); %t wavenumber spectra
ax(10) = axes('Position',[0.52     0.14    0.44      0.3781]); %s wavenumber spectra

box(1).Position = [0.0400    0.0050    0.2262    0.0872];
box(2).Position = [0.2712    0.0050    0.1106    0.0872];
box(3).Position = [0.3869    0.0050    0.1106    0.0872];
box(4).Position = [0.5025    0.0050    0.1106    0.0872];
box(5).Position = [0.6181    0.0050    0.1106    0.0872];
box(6).Position = [0.7337    0.0050    0.2263    0.0872];

%% Define some colors
% ----------------------------------------------------------

cols.t1 = [29 78 140]./255;%[34 63 153]./255;
cols.t2 = [78 173 173]./255;
cols.s1 = [60 134 76]./255;%[235 64 61]./255;
cols.s2 = [173 215 136]./255;%[242 165 59]./255;
cols.a1 = [129 27 112]./255;%[0  0.2667 0.2706];
cols.a2 = [235 64 61]./255;%[0.1725    0.4706    0.4510];
cols.a3 = [245 199 118]./255;%[0.4353    0.7255    0.5608];
cols.T = [185 38 26]./255;
cols.S = [0 0 0];
cols.w = [0 0 0];
cols.profInfo = [174 197 227]./255;
cols.scanInfo = [0.9 0.9 0.9];
cols.batch11 = [0.4902    0.1647    0.4706];
cols.batch12 = [0.6784    0.1529    0.6431];
cols.batch21 = [0.8353    0.1059    0.7922];
cols.batch22 = [0.9679    0.4079    0.6317];
cols.panchev1 = [0.4902    0.1647    0.4706];
cols.panchev2 = [0.8353    0.1059    0.7922];

%     cols.batch11 = [0.9679    0.4079    0.6317];
%     cols.batch12 = [0.8661    0.2040    0.5917];
%     cols.batch21 = [0.6822    0.0051    0.4945];
%     cols.batch22 = [0.4785    0.0039    0.4666];
%     cols.panchev1 = [0.9679    0.4079    0.6317];
%     cols.panchev2 = [0.4785    0.0039    0.4666];

%% Get profile data
% ----------------------------------------------------------

dTdV=[Meta_Data.epsi.t1.dTdV Meta_Data.epsi.t2.dTdV];

% noise floor
logf=log10(Profile.fe);
h_freq=get_filters_MADRE(Meta_Data,Profile.fe);
switch Meta_Data.MAP.temperature
    case 'Tdiff'
        FPO7noise=load(fullfile(Meta_Data.CALIpath,'FPO7_noise.mat'),'n0','n1','n2','n3');
    otherwise
        FPO7noise=load(fullfile(Meta_Data.CALIpath,'FPO7_notdiffnoise.mat'),'n0','n1','n2','n3');
end
n0=FPO7noise.n0; n1=FPO7noise.n1; n2=FPO7noise.n2; n3=FPO7noise.n3;
tnoise=10.^(n0+n1.*logf+n2.*logf.^2+n3.*logf.^3);

shearnoise=load(fullfile(Meta_Data.CALIpath,'shear_noise.mat'),'n0s','n1s','n2s','n3s');
n0s=shearnoise.n0s;
n1s=shearnoise.n1s;
n2s=shearnoise.n2s;
n3s=shearnoise.n3s;
snoise=10.^(n0s+n1s.*logf+n2s.*logf.^2+n3s.*logf.^3);

%% Plot profiles and add profile info text boxes
% ----------------------------------------------------------

% ax(1) - Plot chi profiles
axes(ax(1))
p1 = plot(Profile.chi,Profile.pr);
p1(1).Color = cols.t1;
p1(2).Color = cols.t2;
[p1(:).LineWidth] = deal(2);
ax(1).XScale = 'log';
ax(1).YLabel.String = 'Depth (m)';
ax(1).XLabel.String = '\chi (K^2 s^{-1})';

% ax(2) - Plot epsilon profiles
axes(ax(2))
p2 = plot(Profile.epsilon,Profile.pr);
p2(1).Color = cols.s1;
p2(2).Color = cols.s2;
[p2(:).LineWidth] = deal(2);
ax(2).XScale = 'log';
ax(2).XLabel.String = '\epsilon (W kg^{-1})';

% ax(3) - Plot T and S profiles
[ax3,p3(1),p3(2)]=plotxx(Profile.s,Profile.pr,Profile.t,Profile.pr,{'',''},{'',''},ax(3));
[p3(:).LineWidth] = deal(2);
p3(1).Color = cols.S;
p3(2).Color = cols.T;
ax3(1).XColor = cols.S;
ax3(2).XColor = cols.T;
ax3(1).XLabel.String = 'S';
ax3(2).XLabel.String = 'T (°C)';

% ax(4) - Plot fall speed profile
p4 = plot(ax(4),Profile.w,Profile.pr);
p4.Color = cols.w;
p4.LineWidth = 2;
ax(4).XLabel.String = 'w (m s^{-1})';

% Adjust axes properties
[ax(:).YDir] = deal('reverse');
[ax3(:).YDir] = deal('reverse');
[ax([2,4]).YTickLabel] = deal('');
[ax3(:).YTickLabel] = deal('');

%% Add profile info to text boxes
% ----------------------------------------------------------

% Profile info1
lM = length(Meta_Data.path_mission);
idxSlash = strfind(Meta_Data.path_mission,'/');

% There have been a few iterations of the Meta_Data strucutre. MADRE SN and
% rev stored in different places depending on Meta_Data version.
if isfield(Meta_Data,'Hardware') 
   madreSN = Meta_Data.Hardware.SOM.SN;
   madreRev = Meta_Data.Hardware.SOM.rev;
elseif isfield(Meta_Data,'MADRE')
    madreSN = Meta_Data.MADRE.SN;
    madreRev = Meta_Data.MADRE.rev;
end

if lM>61
    s = find(idxSlash<=61,1,'last');
    annotation('textbox',...
        box(1).Position,...
        'String',{strrep([Meta_Data.mission '  -  ' Meta_Data.vehicle_name '  -  ' Meta_Data.deployment],'_','\_'),...
        Meta_Data.path_mission(1:idxSlash(s)),...
        ['  ', Meta_Data.path_mission(idxSlash(s)+1:end)],...
        ['MADRE SN' madreSN ' rev ' madreRev],...
        ['MAP SN'   Meta_Data.MAP.SN   ' rev ' Meta_Data.MAP.rev],...
        [Meta_Data.CTD.name ' ' Meta_Data.CTD.SN],...
        },...
        'FontSize',infoFontSize,...
        'FontName','Monospaced',...
        'LineStyle','-',...
        'EdgeColor','k',...
        'LineWidth',1,...
        'BackgroundColor',cols.profInfo,...
        'Color','k');
elseif lM<=61
    annotation('textbox',...
        box(1).Position,...
        'String',{strrep([Meta_Data.mission '  -  ' Meta_Data.vehicle_name '  -  ' Meta_Data.deployment],'_','\_'),...
        Meta_Data.path_mission(1:61),...
        '',...
        ['MADRE SN' madreSN ' rev ' madreRev],...
        ['MAP SN'   Meta_Data.MAP.SN   ' rev ' Meta_Data.MAP.rev],...
        [Meta_Data.CTD.name ' ' Meta_Data.CTD.SN],...
        },...
        'FontSize',infoFontSize,...
        'FontName','Monospaced',...
        'LineStyle','-',...
        'EdgeColor','k',...
        'LineWidth',1,...
        'BackgroundColor',cols.profInfo,...
        'Color','k');
end

% Profile info2
annotation('textbox',...
    box(2).Position,...
    'String',{
    sprintf('t1 - SN:   %s',Meta_Data.epsi.t1.SN),...
    sprintf('     dTdV: %3.2f ',Meta_Data.epsi.t1.dTdV),...
    sprintf('     %s - %s',Meta_Data.epsi.t1.ADCfilter,Meta_Data.epsi.t1.ADCconf),...
    sprintf('t2 - SN:   %s',Meta_Data.epsi.t2.SN),...
    sprintf('     dTdV: %3.2f',Meta_Data.epsi.t2.dTdV),...
    sprintf('     %s - %s',Meta_Data.epsi.t2.ADCfilter,Meta_Data.epsi.t2.ADCconf),...
    },...
    'FontSize',infoFontSize,...
    'FontName','Monospaced',...
    'LineStyle','-',...
    'EdgeColor','k',...
    'LineWidth',1,...
    'BackgroundColor',cols.profInfo,...
    'Color','k');

%Profile info3
annotation('textbox',...
    box(3).Position,...
    'String',{
    sprintf('s1 - SN: %s', Meta_Data.epsi.s1.SN),...
    sprintf('     Sv: %3.2f', Meta_Data.epsi.s1.Sv),...
    sprintf('     %s - %s', Meta_Data.epsi.s1.ADCfilter,Meta_Data.epsi.s1.ADCconf),...
    sprintf('s2 - SN: %s',Meta_Data.epsi.s2.SN),...
    sprintf('     Sv:%3.2f ',Meta_Data.epsi.s2.Sv),...
    sprintf('     %s - %s', Meta_Data.epsi.s2.ADCfilter,Meta_Data.epsi.s2.ADCconf),...
    },...
    'FontSize',infoFontSize,...
    'FontName','Monospaced',...
    'LineStyle','-',...
    'EdgeColor','k',...
    'LineWidth',1,...
    'BackgroundColor',cols.profInfo,...
    'Color','k');

%% Loop through scans, plot data, and add scan info text box
% ----------------------------------------------------------

dTdV=[Meta_Data.epsi.t1.dTdV,Meta_Data.epsi.t2.dTdV];
Sv=[Meta_Data.epsi.s1.Sv,Meta_Data.epsi.s2.Sv];
Gr=9.81;

for k=scanNum
    
    
    % Get scan data
    scan = get_scan_spectra(Meta_Data,Profile,k);
    
    if isfield(scan,'pr') %if there's data in this scan...
        
        
        %% Get scan data
        % -------------------------------------------
        
        % Make nan arrays for missing scan data
        if ~isfield(scan,'a2')
            scan.a2 = nan(size(scan.a1,1),size(scan.a1,2));
            scan.Pa.a2 = nan(size(scan.Pa.a1,1),size(scan.Pa.a1,2));
            scan.Cu1a.a2 = nan(size(scan.Cu1a.a1,1),size(scan.Cu1a.a1,2));
            scan.Cu2a.a2 = nan(size(scan.Cu2a.a1,1),size(scan.Cu2a.a1,2));
        end
        if ~isfield(scan.chi,'t1')
            scan.chi.t1 = nan(size(scan.chi.t2,1),size(scan.chi.t2,2));
            scan.Ptgk.t1 = nan(size(scan.Ptgk.t2,1),size(scan.Ptgk.t2,2));
            scan.kc.t1 = nan;
        end
        if ~isfield(scan.chi,'t2')
            scan.chi.t2 = nan(size(scan.chi.t1,1),size(scan.chi.t1,2));
            scan.Ptgk.t2 = nan(size(scan.Ptgk.t1,1),size(scan.Ptgk.t1,2));
        end
        
        % Get FPO7 noise
        k_noise=Profile.fe./Profile.w(k);
        noise_t=tnoise.*dTdV(1).^2./h_freq.FPO7(Profile.w(k));
        tnoise_k= (2*pi*k_noise).^2 .* noise_t.*Profile.w(k);        % T1_k spec  as function of k
        
        % Get shear noise
        TFshear=(Sv(1).*Profile.w(k)/(2*Gr)).^2 .* h_freq.shear.* haf_oakey(Profile.fe,Profile.w(k));
        snoise_k= (2*pi*k_noise).^2 .* snoise.*Profile.w(k)./TFshear;        % T1_k spec  as function of k
        
        % Get Batchelor spectrum for each combination of epsilon/chi
        [kbatch11,Pbatch11] = batchelor(scan.epsilon.s1,scan.chi.t1, ...
            scan.kvis,scan.ktemp);
        [kbatch21,Pbatch21] = batchelor(scan.epsilon.s2,scan.chi.t1, ...
            scan.kvis,scan.ktemp);
        
        [kbatch12,Pbatch12] = batchelor(scan.epsilon.s1,scan.chi.t2, ...
            scan.kvis,scan.ktemp);
        [kbatch22,Pbatch22] = batchelor(scan.epsilon.s2,scan.chi.t2, ...
            scan.kvis,scan.ktemp);
        
        % Smooth Tdiff wavenumber spectra
        smTG1 = smoothdata(scan.Ptgk.t1,'movmean',10);
        smTG2 = smoothdata(scan.Ptgk.t2,'movmean',10);
        
        % Smooth shear wavenumber spectra
        smS1 = smoothdata(scan.Psk.s1,'movmean',10);
        smS2 = smoothdata(scan.Psk.s2,'movmean',10);
        
        % Get the coherences between each shear/acceleration pair
        Co11 = scan.Cu1a.a1;
        Co21 = scan.Cu2a.a1;
        Co12 = scan.Cu1a.a2;
        Co22 = scan.Cu2a.a2;
        Co13 = scan.Cu1a.a3;
        Co23 = scan.Cu2a.a3;
        
               
        %% Frequency plots
        % ----------------------------------------------------------
        
        % ax(5) - Velocity spectra (vs frequency)
        axes(ax(5))
        p5(1) = loglog(Profile.fe,scan.P.s1,'color',cols.s1,'linestyle',':');
        hold on
        p5(2) = loglog(Profile.fe,scan.Pv.s1,'color',cols.s1);
        p5(3) = loglog(Profile.fe,scan.P.s2,'color',cols.s2,'linestyle',':');
        p5(4) = loglog(Profile.fe,scan.Pv.s2,'color',cols.s2);
        grid on
        legend({'s1','s1-noncoh','s2','s2-noncoh'},'location','northwest','numcolumns',2);
        ax(5).XLabel.String = '';
        ax(5).YLabel.String = 's^{-2}/Hz';
        
        % ax(7) - Acceleration spectra (vs frequency)
        axes(ax(7))
        p7(1) = loglog(Profile.fe,scan.Pa.a1,'color',cols.a1);
        hold on
        p7(2) = loglog(Profile.fe,scan.Pa.a2,'color',cols.a2);
        p7(3) = loglog(Profile.fe,scan.Pa.a3,'color',cols.a3);
        grid on
        legend([p7(1),p7(2),p7(3)],'a1','a2','a3','location','northwest','numcolumns',2)
        ax(7).XLabel.String = 'Hz';
        ax(7).YLabel.String = 'g^2/Hz';
        
        
        % ax(6) - Plot coherence with shear channel 1 (vs frequency)
        axes(ax(6))
        p6(1) = semilogx(Profile.fe,Co11,'color',cols.a1);
        hold on
        p6(2) = semilogx(Profile.fe,Co12,'color',cols.a2);
        p6(3) = semilogx(Profile.fe,Co13,'color',cols.a3);
        p6(4) = semilogx(Profile.fe,Profile.Cu1a3,'color',cols.s1);
        grid on
        legend('s1a1','s1a2','s1a3','s1a3 full profile','location','northwest','numcolumns',2)
        ax(6).YAxisLocation='right';
        ax(6).YLabel.String = 'Coherence';
        ax(6).XLabel.String = '';
        
        % ax(8) - Plot coherence with shear channel 2 (vs frequency)
        axes(ax(8))
        p8(1) = semilogx(Profile.fe,Co21,'color',cols.a1);
        hold on
        p8(2) = semilogx(Profile.fe,Co22,'color',cols.a2);
        p8(3) = semilogx(Profile.fe,Co23,'color',cols.a3);
        p8(4) = semilogx(Profile.fe,Profile.Cu1a3,'color',cols.s2);
        grid on
        legend('s2a1','s2a2','s2a3','s2a3 full profile','location','northwest','numcolumns',2)
        ax(8).YAxisLocation='right';
        ax(8).YLabel.String = 'Coherence';
        ax(8).XLabel.String = 'Hz';
        
        % Set some common axes properties
        [ax([5,7]).YLim] = deal([1e-10 1e-3]);
        [ax([6,8]).YLim] = deal([0 1]);
        [ax([5,6]).XTickLabel] = deal('');
        [ax(5:8).XLim] = deal(Profile.fe([1 end]));
        [ax(5:8).XTick] = deal([1 10 100]);
        
        %% Wavenumber plots
        % ----------------------------------------------------------
        
        % ax(9) - Plot Tdiff spectra (vs wavenumber)
        axes(ax(9))
        p9(1) = loglog(scan.ke,scan.Ptgk.t1,':','color',cols.t1,'linewidth',2);
        hold on
        p9(2) = loglog(scan.ke,smTG1,'color',cols.t1,'linewidth',3);
        p9(3) = loglog(scan.ke,scan.Ptgk.t2,':','color',cols.t2);
        p9(4) = loglog(scan.ke,smTG2,'color',cols.t2);
        
        % Add Batchelor spectra
        p9(5) = loglog(kbatch11,Pbatch11,'Color',cols.batch11);
        p9(6) = loglog(kbatch12,Pbatch12,'Color',cols.batch12);
        p9(7) = loglog(kbatch21,Pbatch21,'Color',cols.batch21);
        p9(8) = loglog(kbatch22,Pbatch22,'Color',cols.batch22);
        
        % Add noise
        p9(9) = loglog(k_noise,tnoise_k,'k:');
        
        % Add kc
        indkc=find(scan.ke>scan.kc.t1,1,'first');
        p9(10) = scatter(scan.ke(indkc),smTG1(indkc),'filled','d','sizedata',300,'MarkerEdgeColor','k','markerfacecolor',cols.t1,'linewidth',2);
        indkc=find(scan.ke>scan.kc.t2,1,'first');
        p9(11) = scatter(scan.ke(indkc),smTG2(indkc),'filled','p','sizedata',450,'MarkerEdgeColor','k','markerfacecolor',cols.t2,'linewidth',2);
        
        legend('t1','t1smooth','t2','t2smooth',...
            'batch11','batch12','batch21','batch22',...
            'T-noise','t1_{cutoff}','t2_{cutoff}',...
            'location','southwest','numcolumns',3);
        xlim([6e-1 400])
        ylim([1e-10 1e-1])
        grid on
        xlabel('k (cpm)')
        ylabel('\phi^2_{TG} (C^2 m^{-2} / cpm)')
        
        % ax(10) - Plot shear spectra (vs wavenumber)
        axes(ax(10))
        p10(1) = loglog(scan.ke,scan.Psk.s1,':','color',cols.s1);
        hold on
        p10(2) = loglog(scan.ke,smS1,'color',cols.s1);
        p10(3) = loglog(scan.ke,scan.Psk.s2,':','color',cols.s2);
        p10(4) = loglog(scan.ke,smS2,'color',cols.s2);
        p10(5) = loglog(k_noise,snoise_k,'c','linewidth',2);
        % Add kc
        indkc=find(scan.ke>scan.kc.s1,1,'first');
        p10(6) = scatter(scan.ke(indkc),smS1(indkc),'filled','d','sizedata',300,'MarkerEdgeColor','k','markerfacecolor',cols.s1,'linewidth',2);
        indkc=find(scan.ke>scan.kc.s2,1,'first');
        p10(7) = scatter(scan.ke(indkc),smS2(indkc),'filled','p','sizedata',450,'MarkerEdgeColor','k','markerfacecolor',cols.s2,'linewidth',2);
        
        % Add Panchev
        p10(8) = loglog(scan.kpan.s1,scan.Ppan.s1,'Color',cols.panchev1);
        p10(9) = loglog(scan.kpan.s2,scan.Ppan.s2,'Color',cols.panchev2);
        hold on
        legend('s1','s1smooth','s2','s2smooth','noise','s1_{cutoff}','s2_{cutoff}','Panchev1','Panchev2','location','southwest','numcolumns',2);
        xlim([6e-1 400])
        ylim([1e-10 1e-1])
        grid on
        xlabel('k (cpm)')
        ylabel('\phi^2_{shear} (s^{-2} / cpm)')
        
        
        %% ax(1:4) - Shade location of scan on profile plots
        
        % First, reset font size because it affects the x-limits
        drawnow
        [ax(:).FontSize] = deal(axFontSize);
        [ax3(:).FontSize] = deal(axFontSize);
        
        prLongArray = interp1(Profile.ctdtime,Profile.P,Profile.epsitime);
        prInScan = prLongArray(scan.ind_scan);
        y = [nanmin(prInScan) nanmax(prInScan)];
        
        axes(ax(1))
        hold on
        x = [ax(1).XLim(1) ax(1).XLim(2)];
        f(1) = fill(x([1 1 2 2 1]),y([1 2 2 1 1]),'k');
        f(1).FaceAlpha = 0.2;
        f(1).EdgeColor = 'none';
        
        axes(ax(2))
        hold on
        x = [ax(2).XLim(1) ax(2).XLim(2)];
        f(2) = fill(x([1 1 2 2 1]),y([1 2 2 1 1]),'k');
        f(2).FaceAlpha = 0.2;
        f(2).EdgeColor = 'none';
        
        axes(ax3(2))
        hold on
        x = [ax3(2).XLim(1) ax3(2).XLim(2)];
        f(3) = fill(x([1 1 2 2 1]),y([1 2 2 1 1]),'k');
        f(3).FaceAlpha = 0.2;
        f(3).EdgeColor = 'none';
        
        axes(ax(4))
        hold on
        x = [ax(4).XLim(1) ax(4).XLim(2)];
        f(4) = fill(x([1 1 2 2 1]),y([1 2 2 1 1]),'k');
        f(4).FaceAlpha = 0.2;
        f(4).EdgeColor = 'none';
        
        
        % Scan info1
        annotation('textbox',...
            box(4).Position,...
            'String',{datestr(Profile.dnum(k)),...
            strrep([Meta_Data.mission '  -  ' Meta_Data.vehicle_name '  -  ' Meta_Data.deployment],'_','\_'),...
            sprintf('profile %03.0f',profNum),...
            sprintf('scan %03.0f',k),...
            },...
            'FontSize',infoFontSize,...
            'FontName','Monospaced',...
            'LineStyle','-',...
            'EdgeColor','k',...
            'LineWidth',1,...
            'BackgroundColor',cols.scanInfo,...
            'Color','k');
        
        % Scan info2
        annotation('textbox',...
            box(5).Position,...
            'String',{sprintf('pressure    = %3.1f db',Profile.pr(k)),...
            sprintf('speed       = %1.2f m/s',Profile.w(k)),...
            sprintf('temperature = %2.2f °C',Profile.t(k)),...
            sprintf('salinity    = %2.2f psu',Profile.s(k)),...
            },...
            'FontSize',infoFontSize,...
            'FontName','Monospaced',...
            'LineStyle','-',...
            'EdgeColor','k',...
            'LineWidth',1,...
            'BackgroundColor',cols.scanInfo,...
            'Color','k');
        
        % Scan info3
        annotation('textbox',...
            box(6).Position,...
            'String',{
            sprintf('kinematic viscosity = %1.1e m^2 s^{-1}',scan.kvis),...
            sprintf('scalar diffusivity  = %1.1e m^2 s^{-1}',scan.ktemp),...
            ['\epsilon_{1,2}' sprintf(' = %1.2e, %1.2e (W kg^{-1})',scan.epsilon.s1,scan.epsilon.s2)],...
            ['\chi_{1,2}'     sprintf(' = %1.2e, %1.2e (K^2 s^{-1})',scan.chi.t1,scan.chi.t2)],...
            },...
            'FontSize',infoFontSize,...
            'FontName','Monospaced',...
            'LineStyle','-',...
            'EdgeColor','k',...
            'LineWidth',1,...
            'BackgroundColor',cols.scanInfo,...
            'Color','k');
        
        
        %% Adjust axes properties
        % ----------------------------------------------------------
        drawnow
        
        
        ax(1).XLabel.Units = 'normalized';
        ax(1).XLabel.Position(2) = -0.05;
        ax(2).XLabel.Units = 'normalized';
        ax(2).XLabel.Position(2) = -0.05;
        ax3(1).XLabel.Units = 'normalized';
        ax3(1).XLabel.Position(2) = -0.05;
        ax3(2).XLabel.Units = 'normalized';
        ax3(2).XLabel.Position(2) = 0.94;
        ax(4).XLabel.Units = 'normalized';
        ax(4).XLabel.Position(2) = -0.05;
        
        ax(5).YLabel.Units = 'normalized';
        ax(5).YLabel.Position(1) = -0.1;
        ax(7).YLabel.Units = 'normalized';
        ax(7).YLabel.Position(1) = -0.1;
        ax(7).XLabel.Units = 'normalized';
        ax(7).XLabel.Position(2) = -0.15;
        ax(8).XLabel.Units = 'normalized';
        ax(8).XLabel.Position(2) = -0.15;
        
        ax(9).XLabel.Units = 'normalized';
        ax(9).XLabel.Position(2) = -0.05;
        ax(9).YLabel.Units = 'normalized';
        ax(9).YLabel.Position(1) = -0.05;
        
        ax(10).XLabel.Units = 'normalized';
        ax(10).XLabel.Position(2) = -0.05;
        ax(10).YLabel.Units = 'normalized';
        ax(10).YLabel.Position(1) = -0.05;
        
        
    else
        p5 = [];
        p6 = [];
        p7 = [];
        p8 = [];
        p9 = [];
        p10 = [];
        f = [];
    end %end if there's data in this scan
    
end %end loop through scans