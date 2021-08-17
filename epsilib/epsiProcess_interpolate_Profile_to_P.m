function [gridP] = interpolateProfileToP(Profile,P)

varList0 = {'dnum','profNum'};
varList1 = {'w','t'};
varList2 = {'epsilon','epsilon_co','chi'};
for iVar=1:length(varList0) 
    eval(['gridP.',varList0{iVar},' =  nanmean(Profile.(varList0{iVar}));']);
end
gridP.P = P;
for iVar=1:length(varList1)
    eval(['gridP.',varList1{iVar},' = interp1(Profile.pr,Profile.(varList1{iVar}),P);']);
end
for iVar=1:length(varList2) 
    eval(['gridP.',varList2{iVar},'1 = interp1(Profile.pr,Profile.(varList2{iVar})(:,1),P);']);
    eval(['gridP.',varList2{iVar},'2 = interp1(Profile.pr,Profile.(varList2{iVar})(:,2),P);']);
end

% Acceleration - find average power spectral density of verical
% acceleration within a frequency range.
f = Profile.f;
f1 = 10;
f2 = 45;
for iScan=1:size(Profile.Pa_g_f.a1,1)
    [PSDavg(iScan,1)] = avg_psd_in_frange(Profile.Pa_g_f.a1(iScan,:),f,f1,f2);
end
eval(['gridP.a1_avg = interp1(Profile.pr,PSDavg,P);']);

% Use altimeter and pressure to get bottom depth
% Interpolate pressure to altimeter time, and use all times when altimeter < 35
if isfield(Profile,'alttime')
    pres = interp1(Profile.ctdtime,Profile.P,Profile.alttime);
    alt = Profile.dst;
    alt(alt>35) = nan;
    bottom_depth = alt+pres;
    % Take the average of the deepest 20 measurements (20 seconds). This step
    % is an attempt to get rid of any spurious readings that might have
    % occurred further up in the profile
    [~,idxDeep] = sort(pres);
    bottom_depth_mean = nanmean(bottom_depth(idxDeep(end-19:end)));
    eval(['gridP.bottom_depth = bottom_depth_mean;']);
end

% %% Grid to isotherms
% 
% clear P
% 
% T = [4:0.02:8.2].';
% P = [];
% EPS1 = [];
% EPS2 = [];
% EPS1CO = [];
% EPS2CO = [];
% DNUM = [];
% W = [];
% PROFNUM = [];
% 
% for iFile=1:length(fileList)
%     
%    load(fullfile(dataDir,fileList(iFile).name));
%    
%    % Find only non-nan values and sort by temperature
%    notNan = find(~isnan(Profile.t));
%    [~,idxSorted] = sort(Profile.t(notNan));
%    idx = notNan(idxSorted);
%    
%    
%    eps1 = interp1(Profile.t(idx),Profile.epsilon(idx,1),T);
%    eps2 = interp1(Profile.t(idx),Profile.epsilon(idx,2),T);
%    eps1co = interp1(Profile.t(idx),Profile.epsilon_co(idx,1),T);
%    eps2co = interp1(Profile.t(idx),Profile.epsilon_co(idx,2),T);
%    dnum = nanmean(Profile.dnum(idx));
%    w = interp1(Profile.t(idx),Profile.w(idx),T);
%    profNum = Profile.profNum;
%    p = interp1(Profile.t(idx),Profile.pr(idx),T);
%    
%    EPS1 = [EPS1,eps1];
%    EPS2 = [EPS2,eps2];
%    EPS1CO = [EPS1CO,eps1co];
%    EPS2CO = [EPS2CO,eps2co];
%    DNUM = [DNUM,dnum];
%    W = [W,w];
%    PROFNUM = [PROFNUM,profNum];
%    P = [P,p];
%    
%    clear Profile
% end
% 
% gridT = struct('dnum',DNUM,...
%                 'profnum',PROFNUM,...
%                 'pres',P,...
%                 'temp',T,...
%                 'w',W,...
%                 'eps1',EPS1,...
%                 'eps2',EPS2,...
%                 'eps1co',EPS1CO,...
%                 'eps2co',EPS2CO);
% 
% save griddedProfiles -append gridT

% %% Make a plot
% pcolor(gridP.dnum,gridP.pres,log10(gridP.eps2co));
% shading flat
% hold on
% [c,ch] = contour(gridP.dnum,gridP.pres,gridP.temp,'k','levellist',3.6:0.05:6);
% clabel(c,ch)
% ch.LabelSpacing = 300;
% set(gca,'ydir','reverse','ylim',[1200,2000],'clim',[-10,-6])
% colorbar
% %datetick(gca,'x','HH:MM')


