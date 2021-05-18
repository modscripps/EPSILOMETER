function mod_epsilometer_batch_process_v2(Meta_Data,profileIdx)

% mod_epsilometer_batch_process loads the CTD and EFE profiles and run
% mod_epsilometer_turbulence.Depending on the vehicle (Fish /WW)
% it selects the down or up cast.
%
% OUTPUTS
%   - no output to Matlab workspace
%   - Profile###.mat structures are saved in Meta_Data.L1path
%
% INPUTS
%   Meta_Data
%   profileIdx - (optional) list of profile to process (indices of
%   CTDProfiles and EpsiProfiles)
%
% Meta_Data contains all the information required to process epsilon and
% chi (e.g., Path to data ,path to library, length in second of the scans,
% the number and the names of the EFE channels, CTD and EFE sampling
% frequency). Beside the path, a default Meta_Data can be generated with
% mod_epsilometer_create_Meta_Data.m
%
%  Created by Arnaud Le Boyer on 7/28/18.
%  Copyright � 2018 Arnaud Le Boyer. All rights reserved.

%% Download EPSI and CTD profile
load(fullfile(Meta_Data.L1path,['Profiles_' Meta_Data.deployment '.mat']),'CTDProfiles','EpsiProfiles');
if ~exist('CTDProfiles','var')
    load(fullfile(Meta_Data.L1path,['Profiles_' Meta_Data.deployment '.mat']),'CTDProfile','EpsiProfile');
    CTDProfiles = CTDProfile;
    EpsiProfiles = EpsiProfile;
end

switch Meta_Data.vehicle_name
    case 'FISH'
        CTD_Profiles = CTDProfiles.datadown;
        EPSI_Profiles = EpsiProfiles.datadown;
    case {'WW','Seacycler'}
        CTD_Profiles = CTDProfiles.dataup;
        EPSI_Profiles = EpsiProfiles.dataup;
    otherwise
        CTD_Profiles = CTDProfiles.datadown;
        EPSI_Profiles = EpsiProfiles.datadown;
end

%% If no profile is specified, process them all
if nargin<2
    profileIdx = 1:length(CTD_Profiles);
end

%% Parameters fixed by data structure
% length of 1 scan in second
tscan      =   Meta_Data.PROCESS.tscan;
Fs_epsi    =   Meta_Data.AFE.FS;


% add pressure from ctd to the epsi profile. This should be ter mporary until
% the addition of the pressure sensor on Epsi
L = tscan*Fs_epsi;
sav_var_name = [];
nb_profile_perfile = 0;

[~,Meta_Data.PROCESS.fe]  =  pwelch(0*(1:Meta_Data.PROCESS.nfft),...
    Meta_Data.PROCESS.nfft,[], ...
    Meta_Data.PROCESS.nfft, ...
    Meta_Data.AFE.FS,'psd');

Meta_Data.PROCESS.h_freq = get_filters_SOM(Meta_Data,Meta_Data.PROCESS.fe);

% get FPO7 channel average noise to compute chi
switch Meta_Data.AFE.temp_circuit
    case 'Tdiff'
        Meta_Data.PROCESS.FPO7noise = load(fullfile(Meta_Data.CALIpath,'FPO7_noise.mat'),'n0','n1','n2','n3');
    otherwise
        Meta_Data.PROCESS.FPO7noise = load(fullfile(Meta_Data.CALIpath,'FPO7_notdiffnoise.mat'),'n0','n1','n2','n3');
end


for i = profileIdx
    fprintf('Profile %i of %i, (batch processing profiles %i - %i)\n',i,length(EPSI_Profiles),profileIdx(1),profileIdx(end));

    % some cleaning on the Pressure channel
    % w = dPdt needs to be smooth.
    CTD_Profiles{i}.P = filloutliers(CTD_Profiles{i}.P,'center','movmedian',1000);

    CTD_Profiles{i} = structfun(@(x) fillmissing(x,'linear'),CTD_Profiles{i},'Un',0);
    EPSI_Profiles{i} = structfun(@(x) fillmissing(double(x),'linear'),EPSI_Profiles{i},'Un',0);

    %in case there is a mismatch between ctd and epsi time which still
    %happens as of April 17th 2020.
    EPSI_Profiles{i}.epsitime = EPSI_Profiles{i}.epsitime+ ...
        (CTD_Profiles{i}.ctdtime(1)-EPSI_Profiles{i}.epsitime(1));

    %% TODO check the dimension of the raw time series to remove this line
    % profile to process
    Prmin = min(CTD_Profiles{i}.P);
    Prmax = max(CTD_Profiles{i}.P);
    Profile = mod_epsilometer_merge_profile(Meta_Data,CTD_Profiles{i},EPSI_Profiles{i},Prmin,Prmax);
    
    Profile = mod_epsilometer_calc_turbulence_v2(Meta_Data,Profile);
    
    % Add profNum to the top of the Profile structure
    fieldOrder = fields(Profile);
    Profile.profNum = i;
    Profile = orderfields(Profile,['profNum';fieldOrder]);
    
    % Save Profile
    save_var_name = 'Profile';
    save_file_name = sprintf('Profile%03i',i);
    save_file = fullfile(Meta_Data.L1path, ...
                        [save_file_name '.mat']);
    eval(['save(''' save_file ''', ''' save_var_name ''');']);

end
