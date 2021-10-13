function [grid] = epsiProcess_interpolate_Profile_to_P(Profile,P)

grid.mission = Profile.Meta_Data.mission;
grid.vehicle_name = Profile.Meta_Data.vehicle_name;
grid.deployment = Profile.Meta_Data.deployment;

varList0 = {'profNum','dnum','latitude','longitude'};
varList1 = {'w','t','s','th','sgth'};
varList2 = {'epsilon','epsilon_co','chi'};

% Add varList0
for iVar=1:length(varList0) 
    if isfield(Profile,varList0{iVar})
    eval(['grid.',varList0{iVar},' =  nanmean(Profile.(varList0{iVar}));']);
    end
end

% Add depth field and interpolate
grid.pr = P(:);
try
    grid.z = sw_dpth(grid.pr,nanmean(grid.latitude));
catch
    grid.z = sw_dpth(grid.pr,Profile.Meta_Data.PROCESS.latitude);
end

% Add varList1
notNan = ~isnan(Profile.pr);
for iVar=1:length(varList1)
    if isfield(Profile,varList1{iVar})
    grid.(varList1{iVar}) = interp1(Profile.pr(notNan),Profile.(varList1{iVar})(notNan),grid.pr);
    end
end

% Add varList2
for iVar=1:length(varList2) 
    eval(['grid.',varList2{iVar},'1 = interp1(Profile.pr(notNan),Profile.(varList2{iVar})(notNan,1),grid.pr);']);
    eval(['grid.',varList2{iVar},'2 = interp1(Profile.pr(notNan),Profile.(varList2{iVar})(notNan,2),grid.pr);']);
end

% Acceleration - find average power spectral density of verical
% acceleration within a frequency range.
f = Profile.f;
f1 = 10;
f2 = 45;
for iScan=1:size(Profile.Pa_g_f.a1,1)
    [PSDavg(iScan,1)] = avg_psd_in_frange(Profile.Pa_g_f.a1(iScan,:),f,f1,f2);
end
grid.a1_avg = interp1(Profile.pr(notNan),PSDavg(notNan),grid.pr);

% Use altimeter and pressure to get bottom depth
% Interpolate pressure to altimeter time, and use all times when altimeter < 35
if isfield(Profile,'alt')
    [~,iU] = unique(Profile.ctd.time_s);
    ctdZ = interp1(Profile.ctd.time_s(iU),Profile.ctd.z(iU),Profile.alt.time_s);
    hab = Profile.alt.dst;
    hab(hab>35) = nan;
    bottom_depth = hab+ctdZ;
    % Take the average of the deepest 20 measurements (20 seconds). This step
    % is an attempt to get rid of any spurious readings that might have
    % occurred further up in the profile
    [~,idxDeep] = sort(ctdZ);
    bottom_depth_mean = nanmean(bottom_depth(idxDeep(end-19:end)));
    grid.bottom_depth = bottom_depth_mean;
end

