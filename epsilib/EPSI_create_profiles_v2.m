function EPSI_create_profiles_v2(Meta_Data,min_depth,crit_depth,ask_for_input)
% New version of EPSI_create_profiles, modified to work with epsi_class.m
%
%function EPSI_create_profiles(Meta_Data,crit_speed,crit_filt)
% July 2018 ALB
%  split times series into profiles
%
%  input:
% . Meta_Data
%  created with Meta_Data=create_Meta_Data(file). Meta_Data contain the
%  path to calibration file and EPSI configuration needed to process the
%  epsi data
%   
% . crit_speed: speed criterium to to define starts and ends of a cast
% . crit_filt: number in sample for filter the fall rate (dP/dt) time serie
%  this use to filter out the small fluctuation in the fall rate. 
%  It should be changed when the deployment is not regular.  
%  crit speed and crit_filt are used in EPSI_getcastCTD.
%
%  ask_for_input (optional): 0=no input, save profiles automatically
%                            1=ask for input before saving profiles
%  Created by Arnaud Le Boyer on 7/28/18.

%
if nargin==1
    min_depth=3;
    crit_depth=10;
    ask_for_input=1;
end

if nargin==3
    ask_for_input=1;
end

% load CTD and EPSI data
if exist(fullfile(Meta_Data.CTDpath,['ctd_' Meta_Data.deployment '.mat']),'file')
%     CTD=load(fullfile(Meta_Data.CTDpath,['ctd_' Meta_Data.deployment '.mat']),['ctd_' Meta_Data.deployment]'aux1time','T','P','S','sig');
    CTD=load(fullfile(Meta_Data.CTDpath,['ctd_' Meta_Data.deployment '.mat']),'ctd');
    CTD=CTD.ctd;
    % remove nans from the raw CTD data
%     indOK=~isnan(CTD.aux1time);
%     CTD=structfun(@(x) x(indOK),CTD,'un',0);
%NC added the following because SODA WW d3 CTD structure had aux1time but no ctdtime
if any(strcmp(fields(CTD),'aux1time')) && ~any(strcmp(fields(CTD),'ctdtime'))
    CTD.ctdtime=CTD.aux1time;
end
    
    % define casts using CTD data
    [CTDProfiles.up,CTDProfiles.down,CTDProfiles.dataup,CTDProfiles.datadown] = ...
        mod_getcastctd_v2(Meta_Data,min_depth,crit_depth);
else
    CTD=load(fullfile(Meta_Data.CTDpath,['ctd_' Meta_Data.deployment '.mat']),'ctd');
    CTD=CTD.ctd;
    load(fullfile(Meta_Data.L1path,['Profiles_' Meta_Data.deployment '.mat']),'CTDProfiles')
end
EPSI=load(fullfile(Meta_Data.Epsipath,['epsi_' Meta_Data.deployment '.mat']));
EPSI = EPSI.epsi;


%% plot pressure and highlight up /down casts in red/green
if min(CTD.P)>8
    warning('there is an offset on the RBR pressure. We substract 8 m to get the profile.');
    CTD.P=CTD.P-min(CTD.P);
end

close all
plot(CTD.ctdtime,CTD.P)
hold on
for i=1:length(CTDProfiles.up)
       plot(CTDProfiles.dataup{i}.ctdtime,CTDProfiles.dataup{i}.P,'r')
end
for i=1:length(CTDProfiles.down)
       plot(CTDProfiles.datadown{i}.ctdtime,CTDProfiles.datadown{i}.P,'g')
end

set(gca,'xlim',[nanmin(CTDProfiles.datadown{1}.ctdtime)-0.01,nanmax(CTDProfiles.dataup{end}.ctdtime)+0.01])
datetick('x','keeplimits')

%MHA: plot ocean style
axis ij
figureStamp(mfilename('fullpath'))
print('-dpng2',[Meta_Data.CTDpath 'Profiles_Pr.png'])


%% do we want to save or change the speed and filter criteria
str1 = strvcat(sprintf('--- User-defined minimum profile depth is %2.1f m and minimum profile length is %3.1f',min_depth,crit_depth));
str2 = strvcat('--- User did not define minimum profile depth or minimum profile length.',...
    sprintf('Default values were used: minimum profile depth is %2.1f m and minimum profile length is %3.1f',min_depth,crit_depth));

if nargin==1
    disp(str2)
elseif nargin>1
    disp(str1)
end

switch ask_for_input
    case 1
        answer1  = input('??? Do you want to keep these parameters and save these profiles? (yes/no)','s');
    case 0
        answer1 = 'yes';
end

switch answer1
    case {'yes','y'}
        [EpsiProfiles.up,EpsiProfiles.down,EpsiProfiles.dataup,EpsiProfiles.datadown] =...
            mod_getcastepsi(EPSI,CTD.ctdtime,CTDProfiles.up,CTDProfiles.down);
                
        filepath=fullfile(Meta_Data.L1path,['Profiles_' Meta_Data.deployment '.mat']);
        fprintf('Saving data in %s \n',filepath)
        save(filepath,'CTDProfiles','EpsiProfiles','-v7.3');
        
        Meta_Data.nbprofileup=numel(CTDProfiles.up);
        Meta_Data.nbprofiledown=numel(CTDProfiles.down);
        Meta_Data.maxdepth=max(cellfun(@(x) max(x.P),CTDProfiles.dataup));
        save(fullfile(Meta_Data.L1path,'Meta_Data.mat'),'Meta_Data')
    case {'no','n'}
        disp(' ')
        disp('*********************')
        disp('User should adjust values for Meta_Data.PROCESS.Prmin_prof and Meta_Data.PROCESS.Prcrit_prof');
        disp(' ')
end




