function mod_epsilometer_grid_plot(Map,Meta_Data)

% plot depth time maps of an epsilometer deployment
% plots epsilon1, epsilon2, chi1, chi2, epsilon from chi1 and epsilon from
% chi2.
% save the plots in L1 folder
% TODO add coherence plots, qc flag, t, s, w, acceleration
% TODO in this plot all axes are on individual figures
% the routine could output the axes so we could change and arrange them
% at will a posteriori

% aleboyer@ucsd.edu 03/04/2020

Map.sgth=filloutliers(sw_dens(Map.s,Map.t,Map.z).','nearest','movmedian',10).';
Map.level_sig=linspace(min(nanmean(Map.sgth,2)),max(nanmean(Map.sgth,2)),100);
Map.eta=zeros(100,numel(Map.dnum));

for dt=1:numel(Map.dnum)
    indnan=~isnan(Map.sgth(:,dt));
    if sum(indnan)>1
    Map.eta(:,dt)=interp1(Map.sgth(indnan,dt),Map.z(indnan),Map.level_sig);
    end
end
dvals2=floor(nanmean(Map.eta,2)./2);
dmeta2=diff(dvals2);
Map.eta2m=Map.eta(dmeta2>0,:);


%eps_chi = chi *N^2 / gamma / T_z^2 where gamma = 0.2

[T,Z]=size(Map.t);
Map.gamma=.2;
zaxis2D=repmat(Map.z,[1,Z]);

% despite tyhe fact that sw_bfrq claims the results is in s^{-2} 
% it is in fact in (rad/s^{-1})^2
%N2 = sw_bfrq(s,t,zaxis2D,[])./(2*pi)^2; 
Map.N2 = sw_bfrq(Map.s,Map.t,zaxis2D,[]); 
%%
Tz=diff(Map.t)./diff(zaxis2D);
%%
zaxis12=Map.z(1:end-1)+diff(Map.z);
chi12=interp1(Map.z,Map.chi1,zaxis12);
chi22=interp1(Map.z,Map.chi2,zaxis12);

Map.epsi_chi1 = interp1(zaxis12,chi12.* Map.N2 ./Map.gamma ./ Tz.^2,Map.z);
Map.epsi_chi2 = interp1(zaxis12,chi22.* Map.N2 ./Map.gamma ./ Tz.^2,Map.z);
Map.epsi_chi1(Map.epsi_chi1<0)=nan;
Map.epsi_chi2(Map.epsi_chi2<0)=nan;

Map.N2=interp1(zaxis12,Map.N2,Map.z);
Map.N2(Map.N2<=0)=nan;
Map.N2=fillmissing(Map.N2,'linear');




% %%
close all

% epsilon 1 
fontsize=25;
figure;
colormap('parula')
pcolor(Map.dnum,Map.z,log10(Map.epsilon_co1));shading flat;axis ij
hold on
plot(Map.dnum,Map.eta2m,'Color',[.1,.1,.1,.6],'linewidth',1)
colorbar
caxis([-10,-5])
set(gca,'XTickLabelRotation',25)
datetick
cax=colorbar;
xlabel(['Start date :' datestr(Map.dnum(1),'mm-dd-yyyy')],'fontsize',fontsize)
set(gca,'fontsize',fontsize)
ylabel(cax,'log_{10}(\epsilon)','fontsize',fontsize)
ylabel('Depth (m)','fontsize',fontsize)

fig=gcf;
fig.PaperPosition = [0 0 12 8];
fig.PaperOrientation='Portrait';
print(fullfile(Meta_Data.paths.profiles,'Epsi1_map.png'),'-dpng2')

% epsilon 2
figure;
colormap('parula')
pcolor(Map.dnum,Map.z,log10(Map.epsilon_co2));shading flat;axis ij
hold on
plot(Map.dnum,Map.eta2m,'k')
colorbar
caxis([-11,-5])
set(gca,'XTickLabelRotation',45)
datetick
cax=colorbar;
xlabel(['Start date :' datestr(Map.dnum(1),'mm-dd-yyyy')],'fontsize',15)
set(gca,'fontsize',15)
ylabel(cax,'log_{10}(\epsilon)','fontsize',20)
ylabel('Depth (m)','fontsize',20)

fig=gcf;
fig.PaperPosition = [0 0 15 10];
fig.PaperOrientation='Portrait';
print(fullfile(Meta_Data.paths.profiles,'Epsi2_Map.png'),'-dpng2')


% chi 1 
figure;
colormap('parula')
pcolor(Map.dnum,Map.z,log10(Map.chi1));shading flat;axis ij
hold on
plot(Map.dnum,Map.eta2m,'Color',[.1,.1,.1,.6],'linewidth',1)
colorbar
caxis([-11,-5])
set(gca,'XTickLabelRotation',25)
datetick
cax=colorbar;
xlabel(['Start date :' datestr(Map.dnum(1),'mm-dd-yyyy')],'fontsize',15)
set(gca,'fontsize',fontsize)
ylabel(cax,'log_{10}(\chi)','fontsize',fontsize)
ylabel('Depth (m)','fontsize',fontsize)

fig=gcf;
fig.PaperPosition = [0 0 12 8];
fig.PaperOrientation='Portrait';
print(fullfile(Meta_Data.paths.profiles,'Chi1_map1.png'),'-dpng2')

%chi2 
figure;
colormap('parula')
pcolor(Map.dnum,Map.z,log10(Map.chi2));shading flat;axis ij
hold on
plot(Map.dnum,Map.eta2m,'Color',[.1,.1,.1,.6],'linewidth',1)
colorbar
caxis([-11,-5])
set(gca,'XTickLabelRotation',45)
datetick
cax=colorbar;
xlabel(['Start date :' datestr(Map.dnum(1),'mm-dd-yyyy')],'fontsize',15)
set(gca,'fontsize',15)
ylabel(cax,'log_{10}(\chi)','fontsize',20)
ylabel('Depth (m)','fontsize',20)

fig=gcf;
fig.PaperPosition = [0 0 15 10];
fig.PaperOrientation='Portrait';
print(fullfile(Meta_Data.paths.profiles,'Chi2_map.png'),'-dpng2')


%epsilon from chi1 
figure;
colormap('parula')
pcolor(Map.dnum,Map.z,log10(Map.epsi_chi1));shading flat;axis ij
hold on
plot(Map.dnum,Map.eta2m,'Color',[.1,.1,.1,.6],'linewidth',1)
colorbar
caxis([-11,-5])
set(gca,'XTickLabelRotation',25)
datetick
cax=colorbar;
xlabel(['Start date :' datestr(Map.dnum(1),'mm-dd-yyyy')],'fontsize',fontsize)
set(gca,'fontsize',fontsize)
ylabel(cax,'log_{10}(\epsilon_{\chi_1})','fontsize',fontsize)
ylabel('Depth (m)','fontsize',fontsize)

fig=gcf;
fig.PaperPosition = [0 0 14 9];
fig.PaperOrientation='Portrait';
print(fullfile(Meta_Data.paths.profiles,'Epsichi1_map.png'),'-dpng2')

%epsilon from chi2 
figure;
colormap('parula')
pcolor(Map.dnum,Map.z,log10(Map.epsi_chi2));shading flat;axis ij
hold on
plot(Map.dnum,Map.eta2m,'Color',[.1,.1,.1,.6],'linewidth',1)
colorbar
caxis([-11,-5])
set(gca,'XTickLabelRotation',45)
datetick
cax=colorbar;
xlabel(['Start date :' datestr(Map.dnum(1),'mm-dd-yyyy')],'fontsize',15)
set(gca,'fontsize',15)
ylabel(cax,'log_{10}(\epsilon_{\chi_2})','fontsize',20)
ylabel('Depth (m)','fontsize',20)

fig=gcf;
fig.PaperPosition = [0 0 15 10];
fig.PaperOrientation='Portrait';
print(fullfile(Meta_Data.paths.profiles,'Epsichi2_map.png'),'-dpng2')

