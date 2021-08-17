% process TEMPLATE mission
%
% -------------------------------------------------------------------------
%% Define process directory and data directory
% process_dir = directory containing epsi library
% configFile  = path to config file
process_dir = '/Volumes/Berry/EPSILOMETER/';
configFile = '/Volumes/Berry/blt_epsi/0710_test/bench_config';
data_dir = fileparts(configFile);
cd(data_dir)
addpath(genpath(process_dir))

%% Initialize epsi class
ec = epsi_class;

%% Pick out profile drops from CTD pressure timeseries
% The seabird on the mooring outputs everything as profiles already, so  we
% need to grab the tMin and tMax from this output.

% Is the CTD sampling frequnecy 16 Hz? No, it's 1 Hz.
ec.Meta_Data.PROCESS.Fs_ctd = 1;

% CTD data directory:
presFileDir = '/Volumes/Berry/blt_epsi/0705_process_mooring/mp1_ctd_for_epsi';

fileList = dir(fullfile(presFileDir,'mp*.mat'));
for iF=1:length(fileList)
    disp(num2str(iF))
    ctd = load(fullfile(presFileDir,fileList(iF).name));
    tMin = ctd.dnum(1);
    tMax = ctd.dnum(end);

    try
        Profile = ec.f_cropTimeseries(tMin,tMax);
        
        % Add CTD data to Profile
        Profile.ctddnum = ctd.dnum(:);
        Profile.ctdtime = ctd.dnum(:).*24*60*60;
        Profile.P = ctd.p(:);
        Profile.T = ctd.t(:);
        Profile.dPdt = [0; diff(Profile.P)./diff(Profile.ctdtime)];
        
        % Add profile number
        Profile.profNum = iF;
        
        saveName = fullfile(ec.Meta_Data.L1path,sprintf('Profile%03.0f',iF));
        eval(['save ' saveName ' Profile'])
        clear Profile
    catch
        disp(['no profile data for ' num2str(iF)])
    end

end


%% Calibrate FPO7 temperature
% This step finds the longest profile in the deployment and uses it to
% calibrate FPO7 temperature to Seabird temperature. The calibration values
% are saved in Meta_Data.AFE.t1 and .t2
%
% *** ec.f_calibrateTemperature looks for the longest profile in the
% PressureTimeseries strucutres and grabs the appropriate .mat files that
% contain epsi and ctd strucutres to do the calibration. Here, with the ctd
% data separate in the mooring deployment, we'll just pick out one of our
% manually merged profiles to do the calibration on.
load(fullfile(ec.Meta_Data.MATpath,'Epsi_PressureTimeseries'))
profLengths = PressureTimeseries.endprof-PressureTimeseries.startprof;
pRange = PressureTimeseries.P(PressureTimeseries.endprof) -...
    PressureTimeseries.P(PressureTimeseries.startprof);
[~,idxProf] = max(pRange);
load(fullfile(ec.Meta_Data.L1path,sprintf('Profile%03.0f',idxProf)));

ec.Meta_Data = mod_epsi_temperature_spectra_v3(ec.Meta_Data,Profile);

% Check calibration values
ec.Meta_Data.AFE.t1
ec.Meta_Data.AFE.t2

%% Compute turbulence variables
for iF=1:length(fileList)
    load(fullfile(ec.Meta_Data.L1path,sprintf('Profile%03.0f',iF)));
    Profile = ec.f_computeTurbulence(Profile);
    clear Profile
end

% %% Compute turbulence variables and grid to P array (optional)
% gridData = 0;
% P = 1000:1:2200;
% 
% % Interpolate this profile to standard pressure grid
% gridPnew = obj.f_interpolateProfileToP(Profile,P);
% 
% % Try loading griddedProfiles. If it doesn't exist, we'll
% % make it
% if exist(fullfile(obj.Meta_Data.L1path,'griddedProfiles'),'file')
%     load(fullfile(obj.Meta_Data.L1path,'griddedProfiles'));
% end
% 
% % Add gridded profile to full grid. If full grid doesn't
% % exist yet, create it.
% if exist('gridP','var')
%     varList = fields(gridP);
%     for iVar=1:length(varList)
%         gridP.(varList{iVar})(:,end+1) = gridPnew.(varList{iVar});
%         
%     end
% else
%     varList = fields(gridPnew);
%     for iVar=1:length(varList)
%         gridP.(varList{iVar})(:,1) = gridPnew.(varList{iVar});
%     end
% end
% saveName = fullfile(obj.Meta_Data.L1path,'griddedProfiles');
% eval(['save ' saveName ' gridP']);








