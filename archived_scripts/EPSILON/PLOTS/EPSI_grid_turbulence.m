function EPSI_grid_turbulence(Meta_Data,MS)

% *call datenum vector dnum not time
% *include a yday field
% *all matrices are [depth x time] not the other way around
% *z vector is called z and are column vectors not row vectors
% *n2 not N2
% *sgth not Sig
% *s not S
% *t not T
% *lat and lon and H fields the same dimensions as yday (water depth H can be the same value for all)
% *Then yes, include an info structure with metadata and the processing script name and date.



MSempty=cellfun(@isempty,MS);
Map_pr=cellfun(@(x) (x.pr),MS(~MSempty),'un',0);
z=(min([Map_pr{:}]):.5:max([Map_pr{:}])).';
epsilon2=cellfun(@(x) interp1(x.pr,x.epsilon(:,2),z),MS(~MSempty),'un',0);
epsilon1=cellfun(@(x) interp1(x.pr,x.epsilon(:,1),z),MS(~MSempty),'un',0);
chi1=cellfun(@(x) interp1(x.pr,x.chi(:,1),z),MS(~MSempty),'un',0);
chi2=cellfun(@(x) interp1(x.pr,x.chi(:,2),z),MS(~MSempty),'un',0);
t=cellfun(@(x) interp1(x.pr,x.t,z),MS(~MSempty),'un',0);
s=cellfun(@(x) interp1(x.pr,x.s,z),MS(~MSempty),'un',0);
w=cellfun(@(x) interp1(x.pr,x.w,z),MS(~MSempty),'un',0);
dnum=cell2mat(cellfun(@(x) mean(x.time),MS(~MSempty),'un',0));

epsilon1=cell2mat(epsilon1);
epsilon2=cell2mat(epsilon2);
chi1=cell2mat(chi1);
chi2=cell2mat(chi2);
t=cell2mat(t);
s=cell2mat(s);
w=cell2mat(w);


sgth=filloutliers(sw_dens(s,t,z).','nearest','movmedian',10).';

level_sig=linspace(min(nanmean(sgth,2)),max(nanmean(sgth,2)),100);
for dt=1:numel(dnum)
    indnan=~isnan(sgth(:,dt));
    eta(:,dt)=interp1(sgth(indnan,dt),z(indnan),level_sig);
end
dvals2=floor(nanmean(eta,2)./2);
dmeta2=diff(dvals2);
eta2m=eta(dmeta2>0,:);
plot(eta2m.')


if isfield('Meta_Data','lat')
    lat=dnum*0+lat;
else
    lat=dnum*nan;
end

if isfield('Meta_Data','lon')
    lon=dnum*0+lon;
else
    lon=dnum*nan;
end

if isfield('Meta_Data','H')
    H=dnum*0+H;
else
    H=dnum*nan;
end


%eps_chi = chi *N^2 / gamma / T_z^2 where gamma = 0.2

[T,Z]=size(t);
gamma=.2;
zaxis2D=repmat(z,[1,Z]);

% despite tyhe fact that sw_bfrq claims the results is in s^{-2} 
% it is in fact in (rad/s^{-1})^2
%N2 = sw_bfrq(s,t,zaxis2D,[])./(2*pi)^2; 
N2 = sw_bfrq(s,t,zaxis2D,[]); 
%%
Tz=diff(t)./diff(zaxis2D);
%%
zaxis12=z(1:end-1)+diff(z);
chi12=interp1(z,chi1,zaxis12);
chi22=interp1(z,chi2,zaxis12);


epsi_chi1 = interp1(zaxis12,chi12.* N2 ./gamma ./ Tz.^2,z);
epsi_chi2 = interp1(zaxis12,chi22.* N2 ./gamma ./ Tz.^2,z);
epsi_chi1(epsi_chi1<0)=nan;
epsi_chi2(epsi_chi2<0)=nan;

N2=interp1(zaxis12,N2,z);
N2(N2<=0)=nan;
N2=fillmissing(N2,'linear');
n2=N2;




save(fullfile(Meta_Data.L1path,'Turbulence_grid.mat'), ...
    'epsilon1','epsilon2', ...
    'chi1','chi2', ...
    'dnum','z', ...
    'sgth','t','s','w','eta2m','lat','lon','H','epsi_chi1','epsi_chi2','n2')


% %%
close all

% epsilon 1 
fontsize=25;
figure;
colormap('parula')
% pcolor(dnum,z,log10(real(epsilon1)));shading flat;axis ij
pcolor(dnum,z,sh_qcflag1);shading flat;axis ij
hold on
plot(dnum,eta2m,'Color',[.1,.1,.1,.6],'linewidth',1)
colorbar
% caxis([-10,-5])
set(gca,'XTickLabelRotation',25)
datetick
cax=colorbar;
xlabel(['Start date :' datestr(dnum(1),'mm-dd-yyyy')],'fontsize',fontsize)
set(gca,'fontsize',fontsize)
ylabel(cax,'log_{10}(\epsilon)','fontsize',fontsize)
ylabel('Depth (m)','fontsize',fontsize)

fig=gcf;
fig.PaperPosition = [0 0 12 8];
fig.PaperOrientation='Portrait';
print(fullfile(Meta_Data.L1path,'Epsi1_map.png'),'-dpng2')

% epsilon 2
figure;
colormap('parula')
pcolor(dnum,z,log10(real(epsilon2)));shading flat;axis ij
hold on
plot(dnum,eta2m,'k')
colorbar
caxis([-11,-5])
set(gca,'XTickLabelRotation',45)
datetick
cax=colorbar;
xlabel(['Start date :' datestr(dnum(1),'mm-dd-yyyy')],'fontsize',15)
set(gca,'fontsize',15)
ylabel(cax,'log_{10}(\epsilon)','fontsize',20)
ylabel('Depth (m)','fontsize',20)

fig=gcf;
fig.PaperPosition = [0 0 15 10];
fig.PaperOrientation='Portrait';
print(fullfile(Meta_Data.L1path,'Epsi2_Map.png'),'-dpng2')


% chi 1 
figure;
colormap('parula')
pcolor(dnum,z,log10(real(chi1)));shading flat;axis ij
hold on
plot(dnum,eta2m,'Color',[.1,.1,.1,.6],'linewidth',1)
colorbar
caxis([-11,-5])
set(gca,'XTickLabelRotation',25)
datetick
cax=colorbar;
xlabel(['Start date :' datestr(dnum(1),'mm-dd-yyyy')],'fontsize',15)
set(gca,'fontsize',fontsize)
ylabel(cax,'log_{10}(\chi)','fontsize',fontsize)
ylabel('Depth (m)','fontsize',fontsize)

fig=gcf;
fig.PaperPosition = [0 0 12 8];
fig.PaperOrientation='Portrait';
print(fullfile(Meta_Data.L1path,'Chi1_map1.png'),'-dpng2')

%chi2 
figure;
colormap('parula')
pcolor(dnum,z,log10(real(chi2)));shading flat;axis ij
hold on
plot(dnum,eta2m,'Color',[.1,.1,.1,.6],'linewidth',1)
colorbar
caxis([-11,-5])
set(gca,'XTickLabelRotation',45)
datetick
cax=colorbar;
xlabel(['Start date :' datestr(dnum(1),'mm-dd-yyyy')],'fontsize',15)
set(gca,'fontsize',15)
ylabel(cax,'log_{10}(\chi)','fontsize',20)
ylabel('Depth (m)','fontsize',20)

fig=gcf;
fig.PaperPosition = [0 0 15 10];
fig.PaperOrientation='Portrait';
print(fullfile(Meta_Data.L1path,'Chi2_map.png'),'-dpng2')

%epsilon from chi1 
figure;
colormap('parula')
pcolor(dnum,z,log10(real(epsi_chi1)));shading flat;axis ij
hold on
plot(dnum,eta2m,'Color',[.1,.1,.1,.6],'linewidth',1)
colorbar
caxis([-11,-5])
set(gca,'XTickLabelRotation',25)
datetick
cax=colorbar;
xlabel(['Start date :' datestr(dnum(1),'mm-dd-yyyy')],'fontsize',fontsize)
set(gca,'fontsize',fontsize)
ylabel(cax,'log_{10}(\epsilon_{\chi_1})','fontsize',fontsize)
ylabel('Depth (m)','fontsize',fontsize)

fig=gcf;
fig.PaperPosition = [0 0 14 9];
fig.PaperOrientation='Portrait';
print(fullfile(Meta_Data.L1path,'Epsichi1_map.png'),'-dpng2')

%epsilon from chi2 
figure;
colormap('parula')
pcolor(dnum,z,log10(real(epsi_chi2)));shading flat;axis ij
hold on
plot(dnum,eta2m,'Color',[.1,.1,.1,.6],'linewidth',1)
colorbar
caxis([-11,-5])
set(gca,'XTickLabelRotation',45)
datetick
cax=colorbar;
xlabel(['Start date :' datestr(dnum(1),'mm-dd-yyyy')],'fontsize',15)
set(gca,'fontsize',15)
ylabel(cax,'log_{10}(\epsilon_{\chi_2})','fontsize',20)
ylabel('Depth (m)','fontsize',20)

fig=gcf;
fig.PaperPosition = [0 0 15 10];
fig.PaperOrientation='Portrait';
print(fullfile(Meta_Data.L1path,'Epsichi2_map.png'),'-dpng2')

