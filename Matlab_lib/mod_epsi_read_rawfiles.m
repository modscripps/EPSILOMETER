function mod_epsi_read_rawfiles(Meta_Data)
% Reads the raw .bin files from sd_raw and creates .mat files in three directories:
%   1. Saves SDd**.mat in sd_raw 
%   2. Saves epsi_d**.mat in epsi
%   3. Saves ctd_d**.mat in ctd
% Makes plot check_timestamp.png

% get the filenames of the raw data in the  raw or sdraw folder
switch Meta_Data.PROCESS.recording_mode
    case 'STREAMING'
        list_files = dir(fullfile(Meta_Data.RAWpath,'*.epsi'));
        if isempty(list_files)
            list_files = dir(fullfile(Meta_Data.SDRAWpath,'*.epsi'));
        end
    case 'SD'
        % First, determine which file extension you have, .bin, .dat, or
        % .epsi
        files_in_dir = dir(fullfile(Meta_Data.SDRAWpath));
        isBin = any(contains({files_in_dir.name},'.bin'));
        isDat = any(contains({files_in_dir.name},'.dat'));
        isEpsi = any(contains({files_in_dir.name},'.epsi'));
        if isBin
            list_files = dir(fullfile(Meta_Data.SDRAWpath,'*.bin'));
        elseif isDat
            list_files = dir(fullfile(Meta_Data.SDRAWpath,'*.dat'));
        elseif isEpsi
            list_files = dir(fullfile(Meta_Data.SDRAWpath,'*.epsi'));
        end
end
filenames = {list_files.name};
dirnames = {list_files.folder};

% sort the files from (1 12 112 2 3) to (1 2 3 .... 12 .... 112)
% 
filenames=natsortfiles(filenames);
for i = 1:numel(filenames)
    filenames{i} = fullfile(dirnames{i},filenames{i});
end

% actual read.  be very carefull on the Meta_Data Structure.
a = mod_read_epsi_raw(filenames,Meta_Data);
if ~isfield(a.epsi,'c')
    ramp_count=0*a.epsi.s1_count;
else
    ramp_count=a.epsi.c_count;
end


% save the whole time series 
switch Meta_Data.PROCESS.recording_mode
    case 'STREAMING'

        % mod_epsi_sd_buildtime.m is only used for SD data, and it removes
        % some of the raw and count fields from a. Remove them here, too,
        % so that data formats stay consistent.
        STR.epsi.EPSInbsample = a.epsi.EPSInbsample;
        STR.epsi.t1_volt = a.epsi.t1_volt;
        STR.epsi.t2_volt = a.epsi.t2_volt;
        STR.epsi.s1_volt = a.epsi.s1_volt;
        STR.epsi.s2_volt = a.epsi.s2_volt;
        STR.epsi.c_count = a.epsi.c_count;
        STR.epsi.a1_g = a.epsi.a1_g;
        STR.epsi.a2_g = a.epsi.a2_g;
        STR.epsi.a3_g = a.epsi.a3_g;
        STR.epsi.epsitime=a.epsi.time(:);
        
        STR.aux1.T=a.aux1.T;
        STR.aux1.P=a.aux1.P;
        STR.aux1.C=a.aux1.C;
        STR.aux1.S=sw_salt(a.aux1.C*10./sw_c3515,a.aux1.T,a.aux1.P);
        STR.aux1.sig=sw_pden(STR.aux1.S,a.aux1.T,a.aux1.P,0);
        STR.aux1.aux1time=STR.epsi.epsitime;
        STR.aux1.aux1time=a.aux1.time(:);
        
        a.epsi = STR.epsi;
        a.aux1 = STR.aux1;
        
        save(fullfile(Meta_Data.RAWpath,['STR' Meta_Data.deployment '.mat']),'a','-v7.3')
    case 'SD'
        % the computation of S sig is done in sd buildtime
        %TODO this somewhere else maybe a function common to
        % SD and streaming since there are no good reason to have them separate. 
        SD=mod_epsi_sd_buildtime(Meta_Data,a);
        save(fullfile(Meta_Data.SDRAWpath,['SD' Meta_Data.deployment '.mat']),'a','-v7.3')
        a=SD;
end


% 01/30/2019 we are using c as a ramp samp signal or scan count 
a.epsi.ramp_count=ramp_count;
% -----------------------------
% plot some scan count checks
figure('visible','off');
%ax(1)=subplot(311);plot(1:length(a.madre.EpsiStamp),a.madre.TimeStamp,'b')
ax(1)=axes('position',[0.0700    0.7450    0.9000    0.2050]);
plot(1:length(a.madre.EpsiStamp),a.madre.EpsiStamp,'b')

ax(2)=axes('position',[0.0700    0.5200    0.9000    0.2050]);
plot(1:length(a.madre.EpsiStamp)-1,diff(a.madre.EpsiStamp))
ylim([0 200])

ax(3)=axes('position',[0.0700    0.295    0.9000    0.2050]);
plot((1:length(diff(a.epsi.ramp_count)))/160,diff(a.epsi.ramp_count));
ylim([-5 2])

ax(4)=axes('position',[0.0700    0.07    0.9000    0.2050]);
plot(1:length(a.madre.EpsiStamp),a.madre.TimeStamp/86400+datenum(1970,1,1),'b');
hold on
plot((1:length(a.epsi.epsitime))/160,a.epsi.epsitime,'r');

[ax(:).XTickLabel] = deal('');
linkaxes(ax,'x')
figureStamp(mfilename('fullpath'))
print('-dpng2',fullfile(Meta_Data.SDRAWpath, 'check_timestamp.png'))
close all
% -----------------------------


% save CTD in the ctd folder
if isfield(a,'aux1')
    clear F
    F=fieldnames(a.aux1);
    command=[];
    for f=1:length(F)
        wh_F=F{f}; 
         eval(sprintf('%s.%s=a.aux1.%s;',['ctd_' Meta_Data.deployment ],wh_F,wh_F))
    end
    filepath=fullfile(Meta_Data.CTDpath,['ctd_' Meta_Data.deployment '.mat']);
     save(filepath,['ctd_' Meta_Data.deployment]);
end

% save EPSI in the epsi
clear F
F=fieldnames(a.epsi);
command=[];
for f=1:length(F)
    wh_F=F{f};
    eval(sprintf('%s=a.epsi.%s;',wh_F,wh_F))
    command=[command ',' sprintf('''%s''',F{f})];
end
filepath=fullfile(Meta_Data.Epsipath,['epsi_' Meta_Data.deployment '.mat']);
command=sprintf('save(''%s''%s)',filepath,command);
disp(command)
eval(command);

toc
end


