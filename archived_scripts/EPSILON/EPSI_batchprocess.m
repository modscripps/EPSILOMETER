function EPSI_batchprocess(Meta_Data)

%  input: Meta_Data
%  created with Meta_Data=create_Meta_Data. Meta_Data contain the
%  path to calibration file and EPSI configuration needed to process the
%  epsi data
% e.g .  Meta_Data = 
%          mission: 'For_Arnaud_Sept2018'
%     vehicle_name: 'Granite_unclass_deployments'
%       deployment: 'epsifish4_sep09dep02'
%     path_mission: '/Volumes/DataDrive/GRANITE/'
%          vehicle: 'FISH'
%          process: '~/ARNAUD/SCRIPPS/EPSILOMETER/'
%             root: '/Volumes/DataDrive/GRANITE/For_Arnaud_Sept2018/Granite_unclass_deployments/epsifish4_sep09dep02'
%           L1path: '/Volumes/DataDrive/GRANITE/For_Arnaud_Sept2018/Granite_unclass_deployments/epsifish4_sep09dep02/L1/'
%         Epsipath: '/Volumes/DataDrive/GRANITE/For_Arnaud_Sept2018/Granite_unclass_deployments/epsifish4_sep09dep02/epsi/'
%          CTDpath: '/Volumes/DataDrive/GRANITE/For_Arnaud_Sept2018/Granite_unclass_deployments/epsifish4_sep09dep02/ctd/'
%          RAWpath: '/Volumes/DataDrive/GRANITE/For_Arnaud_Sept2018/Granite_unclass_deployments/epsifish4_sep09dep02/raw/'
%        SDRAWpath: '/Volumes/DataDrive/GRANITE/For_Arnaud_Sept2018/Granite_unclass_deployments/epsifish4_sep09dep02/sd_raw/'
%          PROCESS: [1×1 struct]
%            MADRE: [1×1 struct]
%              MAP: [1×1 struct]
%         Firmware: [1×1 struct]
%             aux1: [1×1 struct]
%             epsi: [1×1 struct]
%         CALIpath: '~/ARNAUD/SCRIPPS/EPSILOMETER/CALIBRATION/ELECTRONICS/'
%           SBEcal: [1×1 struct]
%        starttime: 7.3731e+05
% 
%  Created by Arnaud Le Boyer on 7/28/18.
%  Copyright © 2018 Arnaud Le Boyer. All rights reserved.

% check if already computed the turbulence files. becaue it can be very
% long
if nargin<2
    dsp=0;
end

% download EPSI and CTD profile
try
    load(fullfile(Meta_Data.L1path,['Profiles_' Meta_Data.deployment '.mat']),'CTDProfiles','EpsiProfiles');
catch
    load(fullfile(Meta_Data.L1path,['Profiles_' Meta_Data.deployment '.mat']),'CTDProfile','EpsiProfile');
    CTDProfiles=CTDProfile;
    EpsiProfiles=EpsiProfile;
end
switch Meta_Data.vehicle
    case 'FISH'
        CTD_Profiles=CTDProfiles.datadown;
        EPSI_Profiles=EpsiProfiles.datadown;
    case 'WW'
        CTD_Profiles=CTDProfiles.dataup;
        EPSI_Profiles=EpsiProfiles.dataup;
end
%% Parameters fixed by data structure
% length of 1 scan in second
tscan     =  3;
% Hard cut off frequency to compute espilon. It should below any
% pre-known vibration on the vehicle. For the epsifish the SBE water
% pump is known to vibrate at 50Hz.
Fcut_epsilon=45;  %45 Hz

% sample rate channels
FS        = str2double(Meta_Data.Firmware.sampling_frequency(1:3));
% number of samples per scan (1s) in channels
df        = 1/tscan;

f=(df:df:FS/2)'; % frequency vector for spectra
MS = struct([]);

%TODO interp on a regular time grid epsitime fluctuates.
% this can be imporved with a better us clock in the firmware...
% but it will always varies.
% add pressure from ctd to the epsi profile. This should be ter mporary until
% the addition of the pressure sensor on Epsi
count=0;
rem_nan=@(x) (fillmissing(x,'linear'));
sav_var_name=[];
for i=1:length(EPSI_Profiles)
    L=tscan*320;
    if numel(EPSI_Profiles{i}.epsitime)>10*L
    fprintf('Profile %i over %i\n',i,length(EPSI_Profiles));
    Profile=EPSI_Profiles{i};

    Profile.P=rem_nan(interp1(rem_nan(CTD_Profiles{i}.ctdtime),CTD_Profiles{i}.P,EPSI_Profiles{i}.epsitime));
    Profile.P=filloutliers(Profile.P,'center','movmedian',1000);
    Profile.T=rem_nan(interp1(rem_nan(CTD_Profiles{i}.ctdtime),CTD_Profiles{i}.T,EPSI_Profiles{i}.epsitime));
    Profile.S=rem_nan(interp1(rem_nan(CTD_Profiles{i}.ctdtime),CTD_Profiles{i}.S,EPSI_Profiles{i}.epsitime));
    Profile.time=EPSI_Profiles{i}.epsitime;

    Profile=calc_turbulence(Profile,6,3*325,3*325,325,2,48,Meta_Data);
    eval(sprintf('Profile%03i=Profile;',i))
    clear Profile;
    sav_var_name=[sav_var_name sprintf(',''Profile%03i''',i)];
%     MS{mod(i-1,10)+1}=calc_turbulence(EPSI_Profiles{i},tscan,f,Fcut_epsilon,Meta_Data,dsp,i);
    else
        MS{mod(i-1,10)+1}=[];
    end
    if (mod(i,10)==0)
        save_file=fullfile(Meta_Data.L1path,['Turbulence_Profiles' num2str(count) '.mat']);
        cmd=['save(''' save_file '''' sav_var_name ')'];
        eval(cmd);

        count=count+10;
        sav_var_name=[];
    end
    
end
% if exist('MS','var')
%     save(fullfile(Meta_Data.L1path,['Turbulence_Profiles' num2str(count) '.mat']),cmd(2:end),'-v7.3')
% end






