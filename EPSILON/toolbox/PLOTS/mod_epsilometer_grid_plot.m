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


% %%
close all

% epsilon 1 
fontsize=25;
figure;
colormap('parula')
pcolor(Map.dnum,Map.z,log10(Map.epsilon1));shading flat;axis ij
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
pcolor(Map.dnum,Map.z,log10(Map.epsilon2));shading flat;axis ij
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


% chi 1 
figure;
colormap('parula')
pcolor(Map.dnum,Map.z,log10(Map.chi21));shading flat;axis ij
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
print(fullfile(Meta_Data.paths.profiles,'Chi21_map1.png'),'-dpng2')

%chi2 
figure;
colormap('parula')
pcolor(Map.dnum,Map.z,log10(Map.chi22));shading flat;axis ij
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
print(fullfile(Meta_Data.paths.profiles,'Chi22_map.png'),'-dpng2')



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

