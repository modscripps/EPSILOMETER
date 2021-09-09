function [grid] = epsiProcess_interpolate_Profile_to_z(Profile,z)

grid.mission = Profile.Meta_Data.mission;
grid.vehicle_name = Profile.Meta_Data.vehicle_name;
grid.deployment = Profile.Meta_Data.deployment;


varList0 = {'dnum','profNum'};
varListGPS = {'latitude','longitude'};
varList1 = {'w','t','s','th','sgth'};
varList2 = {'epsilon','epsilon_co','chi'};
for iVar=1:length(varList0) 
    eval(['grid.',varList0{iVar},' =  Profile.(varList0{iVar});']);
end
for iVar=1:length(varListGPS) 
    eval(['grid.',varListGPS{iVar},' =  Profile.gps.(varListGPS{iVar});']);
end
% Add depth field and interpolate
grid.z = z;
for iVar=1:length(varList1)
    eval(['grid.',varList1{iVar},' = interp1(Profile.z,Profile.(varList1{iVar}),z);']);
end
for iVar=1:length(varList2) 
    eval(['grid.',varList2{iVar},'1 = interp1(Profile.z,Profile.(varList2{iVar})(:,1),z);']);
    eval(['grid.',varList2{iVar},'2 = interp1(Profile.z,Profile.(varList2{iVar})(:,2),z);']);
end

% Acceleration - find average power spectral density of verical
% acceleration within a frequency range.
f = Profile.f;
f1 = 10;
f2 = 45;
for iScan=1:size(Profile.Pa_g_f.a1,1)
    [PSDavg(iScan,1)] = avg_psd_in_frange(Profile.Pa_g_f.a1(iScan,:),f,f1,f2);
end
eval(['grid.a1_avg = interp1(Profile.pr,PSDavg,P);']);

% Use altimeter and pressure to get bottom depth
% Interpolate pressure to altimeter time, and use all times when altimeter < 35
if isfield(Profile,'alttime')
    ctdZ = interp1(Profile.ctd.time_s,Profile.ctd.z,Profile.alt.time_s);
    alt = Profile.alt.dst;
    alt(alt>35) = nan;
    bottom_depth = alt+ctdZ;
    % Take the average of the deepest 20 measurements (20 seconds). This step
    % is an attempt to get rid of any spurious readings that might have
    % occurred further up in the profile
    [~,idxDeep] = sort(ctdZ);
    bottom_depth_mean = nanmean(bottom_depth(idxDeep(end-19:end)));
    grid.bottom_depth = bottom_depth_mean);
end

