function Map=mod_epsilometer_grid_turbulence(Meta_Data,MS)

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

Fnames=fieldnames(MS{1});
MSempty=cellfun(@isempty,MS);
for f=1:length(Fnames)
    wh_field=Fnames{f};
    switch wh_field
        case 'pr'
            Map_pr=cellfun(@(x) (x.pr),MS(~MSempty),'un',0);
            Map.z=(min([Map_pr{:}]):.25:max([Map_pr{:}])).';
        case 'epsilon'
            Map.epsilon2=cellfun(@(x) interp1(x.pr,x.epsilon(:,2),Map.z),MS(~MSempty),'un',0);
            Map.epsilon1=cellfun(@(x) interp1(x.pr,x.epsilon(:,1),Map.z),MS(~MSempty),'un',0);
            Map.epsilon1=cell2mat(Map.epsilon1);
            Map.epsilon2=cell2mat(Map.epsilon2);
        case 'chi'
            Map.chi1=cellfun(@(x) interp1(x.pr,x.chi(:,1),Map.z),MS(~MSempty),'un',0);
            Map.chi2=cellfun(@(x) interp1(x.pr,x.chi(:,2),Map.z),MS(~MSempty),'un',0);
            Map.chi1=cell2mat(Map.chi1);
            Map.chi2=cell2mat(Map.chi2);
        case 'chi2'
            Map.chi21=cellfun(@(x) interp1(x.pr,x.chi2(:,1),Map.z),MS(~MSempty),'un',0);
            Map.chi22=cellfun(@(x) interp1(x.pr,x.chi2(:,2),Map.z),MS(~MSempty),'un',0);
            Map.chi21=cell2mat(Map.chi21);
            Map.chi22=cell2mat(Map.chi22);
        case 'dnum'
            Map.dnum=cell2mat(cellfun(@(x) nanmean(x.dnum),MS(~MSempty),'un',0));
        case 'w'
            Map.w=cellfun(@(x) interp1(x.pr,x.w,Map.z),MS(~MSempty),'un',0);
            Map.w=cell2mat(Map.w);
        case 't'
            Map.t=cellfun(@(x) interp1(x.pr,x.t,Map.z),MS(~MSempty),'un',0);
            Map.t=cell2mat(Map.t);
        case 's'
            Map.s=cellfun(@(x) interp1(x.pr,x.s,Map.z),MS(~MSempty),'un',0);
            Map.s=cell2mat(Map.s);
        case 'tg_flag'
            Map.tg_flag1=cellfun(@(x) interp1(x.pr,x.tg_flag(:,1),Map.z),MS(~MSempty),'un',0);
            Map.tg_flag2=cellfun(@(x) interp1(x.pr,x.tg_flag(:,2),Map.z),MS(~MSempty),'un',0);
            Map.tg_flag1=cell2mat(Map.tg_flag1);
            Map.tg_flag2=cell2mat(Map.tg_flag2);
        case 'sh_qcflag'
            Map.sh_qcflag1=cellfun(@(x) interp1(x.pr,x.sh_qcflag(:,1),Map.z),MS(~MSempty),'un',0);
            Map.sh_qcflag2=cellfun(@(x) interp1(x.pr,x.sh_qcflag(:,2),Map.z),MS(~MSempty),'un',0);
            Map.sh_qcflag1=cell2mat(Map.sh_qcflag1);
            Map.sh_qcflag2=cell2mat(Map.sh_qcflag2);
        case {'Pt1','Pt2','Ps1','Ps2','Pa1','Pa2','Pa3','CCu1a1','CCu1a2',...
              'CCu1a3','CCu2a1','CCu2a2','CCu2a3'}  
            Map.(wh_field)=cellfun(@(x) interp1(x.pr,x.(wh_field),Map.z),MS(~MSempty),'un',0);
            Map.(wh_field)=cell2mat(Map.(wh_field));
    end
end


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


if isfield('Meta_Data','lat')
    Map.lat=Map.dnum*0+lat;
else
    Map.lat=Map.dnum*nan;
end

if isfield('Meta_Data','lon')
    Map.lon=Map.dnum*0+lon;
else
    Map.lon=Map.dnum*nan;
end

if isfield('Meta_Data','H')
    Map.H=Map.dnum*0+H;
else
    Map.H=Map.dnum*nan;
end


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

save(fullfile(Meta_Data.paths.profiles,'Turbulence_grid.mat'), ...
    'Map')


