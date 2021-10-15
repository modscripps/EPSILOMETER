function Profile = mod_epsilometer_calc_turbulence(Meta_Data,Profile_or_profNum,saveData)
% Profile = mod_epsilometer_calc_turbulence(Meta_Data,Profile_or_profNum)
%  Profile structure for Micro Structure. Inside Profile you ll find
%  temperature spectra in degC Hz^-1
%  Horizontal  velocity spectra in m^2/s^-2 Hz^-1
%  Acceleration/speed spectra in s^-1 Hz^-1
%
%  Created by Arnaud Le Boyer on 7/28/18.
%  Edited summer 2020 - Nicole Couto

if nargin<3
    saveData = 1;
end


%% Get Profile from Profile_or_profNum
if isnumeric(Profile_or_profNum) && ~isstruct(Profile_or_profNum)
    profNum = Profile_or_profNum;
    load(fullfile(Meta_Data.paths.profiles,sprintf('Profile%03.0f',profNum)));
    
elseif isstruct(Profile_or_profNum)
    Profile = Profile_or_profNum;
elseif isclassfield(Profile_or_profNum,'epsi') && isclassfield(Profile_or_profNum,'ctd') && isclassfield(Profile_or_profNum,'Meta_Data')
    Profile.Meta_Data = Profile_or_profNum.Meta_Data;
    Profile.epsi = Profile_or_profNum.epsi;
    Profile.ctd = Profile_or_profNum.ctd;
else
    error('Need epsi, ctd, and Meta_Data to calculate turbulence parameters!');
end
if isfield(Profile_or_profNum,'gps')
    Profile.latitude = nanmean(Profile_or_profNum.gps.latitude);
    Profile.longitude = nanmean(Profile_or_profNum.gps.longitude);
end

%% despiking shear and temp

%fprintf('despiking Profile\r\n')
%ALB Following the SCOR group recommendation
%I am assuming any impact between probes and stuff in the ocean is about
%50 ms = 50e-3 sec * 320Hz = 16 samples.
%I am removing any outliers that are above  3 times the standard deviation of a
%windows that is 10 * 16 samples. Arbitrarily using 10.
% I am saving the percentage of outliers in the Profile structure
movmean_window_width=10*16;
if ~isempty(Profile.epsi)
    for c=1:Profile.Meta_Data.PROCESS.nb_channels
        wh_channel=Profile.Meta_Data.PROCESS.timeseries{c};
        if ~strcmp(wh_channel,'c_count')
            Profile.epsi.([wh_channel '_raw'])=Profile.epsi.(wh_channel);
            Profile.epsi.(wh_channel)=  ...
                filloutliers(Profile.epsi.([wh_channel '_raw']), ...
                'linear','movmean',movmean_window_width);
            Profile.qc.outliers.(wh_channel)= ...
                sum(Profile.epsi.(wh_channel)~= ...
                Profile.epsi.([wh_channel '_raw']))/...
                numel(Profile.epsi.(wh_channel))*100;
        end
    end
    
end

%% get channels
channels = Meta_Data.PROCESS.timeseries;

nfft = Meta_Data.PROCESS.nfft;
nfftc = Meta_Data.PROCESS.nfftc;
try
    Fs_epsi = Meta_Data.AFE.FS;
catch
    Fs_epsi = Meta_Data.PROCESS.Fs_epsi;
end
Fs_ctd = Meta_Data.PROCESS.Fs_ctd;

fpump = Meta_Data.PROCESS.ctd_fc;
tscan = Meta_Data.PROCESS.tscan;
limit_speed = .2; % limit speed 20 cm s^{-1}
dz  =   Meta_Data.PROCESS.dz;

[~,Meta_Data.PROCESS.fe]  =  pwelch(0*(1:Meta_Data.PROCESS.nfft),...
    Meta_Data.PROCESS.nfft,[], ...
    Meta_Data.PROCESS.nfft, ...
    Fs_epsi,'psd');

if isfield(Meta_Data,'MADRE')
    Meta_Data.PROCESS.h_freq = get_filters_MADRE(Meta_Data,Meta_Data.PROCESS.fe);
else
    Meta_Data.PROCESS.h_freq = get_filters_SOM(Meta_Data,Meta_Data.PROCESS.fe);
end

% get FPO7 channel average noise to compute chi
if isfield(Meta_Data,'MAP')
    temp_circuit = Meta_Data.MAP.temperature;
elseif isfield(Meta_Data,'AFE')
    temp_circuit = Meta_Data.AFE.temp_circuit;
end

switch temp_circuit
    case 'Tdiff'
        Meta_Data.PROCESS.FPO7noise = load(fullfile(Meta_Data.paths.calibration,'FPO7_noise.mat'),'n0','n1','n2','n3');
    otherwise
        Meta_Data.PROCESS.FPO7noise = load(fullfile(Meta_Data.paths.calibration,'FPO7_notdiffnoise.mat'),'n0','n1','n2','n3');
end

%% Cut profile to compute coherence
% Find max and min pressure
Prmin = Meta_Data.PROCESS.Prmin(Profile.ctd.P);
Prmax = Meta_Data.PROCESS.Prmax(Profile.ctd.P);

% Find ctdtime values within this pressure range
inRange = Profile.ctd.P>=Prmin & Profile.ctd.P<=Prmax;
ctdtime = Profile.ctd.time_s(inRange);

% Remove epsi values outside the time period defined by ctdtime
keepEpsi = Profile.epsi.time_s>=nanmin(ctdtime) & Profile.epsi.time_s<=nanmax(ctdtime);
for iChan=1:numel(channels)
    wh_channel = channels{iChan};
    if ~strcmp(wh_channel,'c_count')
        Profile_coh.(wh_channel) = Profile.epsi.(wh_channel)(keepEpsi);
    end
end

% compute coherence with a3 over the full profile.
% NC added 10/14/21 - compute coherence with a1 for earlier datasets. Check
% Meta_Data.PROCESS to determine which channel
if isfield(Meta_Data.PROCESS,'coherent_acc_chan')
    channel = [Meta_Data.PROCESS.coherent_acc_chan '_g'];
else
    channel = 'a3_g';
end
[Profile.Cs1a3_full,Profile.Cs2a3_full,...
    ~,~,~] = mod_efe_scan_coherence(Profile_coh,channel,Meta_Data);


%% define a Pressure axis to an which I will compute epsilon and chi.
%  The spectra will be nfft long centered around P(z) +/- tscan/2.
%
Pr = ceil(min(Profile.ctd.P)):dz:floor(max(Profile.ctd.P));
nbscan = length(Pr);

LCTD = length(Profile.ctd.P);% length of profile
% number of samples for a scan. I make sure it is always even
N_epsi = tscan.*Fs_epsi-mod(tscan*Fs_epsi,2);
N_ctd = tscan.*Fs_ctd-mod(tscan*Fs_ctd,2);

% get shear probe calibration coefficient.
if isfield(Meta_Data,'AFE')
    field_name = 'AFE';
elseif isfield(Meta_Data,'epsi')
    field_name = 'epsi';
end
Sv1 = Meta_Data.(field_name).s1.cal;
Sv2 = Meta_Data.(field_name).s2.cal;
% Sensitivity of FPO7 probe, nominal
dTdV(1) = Meta_Data.(field_name).t1.cal; % define in mod_epsi_temperature_spectra
dTdV(2) = Meta_Data.(field_name).t2.cal; % define in mod_epsi_temperature_spectra

f = Meta_Data.PROCESS.fe;
h_freq = Meta_Data.PROCESS.h_freq;
FPO7noise = Meta_Data.PROCESS.FPO7noise;


% start creating the Profile structure
% get nb of scans in the profile
Profile.pr         =  Pr(:);
Profile.nbscan     =  nbscan;
Profile.nfft       =  nfft;
Profile.nfftc      =  nfftc;
Profile.tscan      =  tscan;
Profile.fpump      =  fpump; % arbitrary cut off frequency usually extract from coherence spectra shear/accel
Profile.f          =  f(:).';

% get index of the acceleration channels (useful when the number of channels is not 8)
idxA = contains(Meta_Data.PROCESS.timeseries,'a');
accList = Meta_Data.PROCESS.timeseries(idxA);



%% Get data for each scan and store averages in Profile

% ------------------------------------------------
% Pre-allocate profile space
Profile.dnum = nan(nbscan,1);
Profile.z = nan(nbscan,1);
Profile.t = nan(nbscan,1);
Profile.w = nan(nbscan,1);
Profile.s = nan(nbscan,1);
Profile.th = nan(nbscan,1);
Profile.sgth = nan(nbscan,1);
Profile.ind_range_ctd = nan(nbscan,2);
Profile.ind_range_epsi = nan(nbscan,2);
Profile.fom = nan(nbscan,2);
Profile.epsilon = nan(nbscan,2);
Profile.epsilon_co = nan(nbscan,2);
Profile.sh_fc = nan(nbscan,2);
Profile.chi = nan(nbscan,2);
Profile.tg_fc = nan(nbscan,2);
Profile.flag_tg_fc = nan(nbscan,2);
Profile.kvis = nan(nbscan,1);
Profile.k = nan(nbscan,length(f));
Profile.Ps_volt_f.s1 = nan(nbscan,length(f));
Profile.Ps_volt_f.s2 = nan(nbscan,length(f));
Profile.Ps_shear_k.s1 = nan(nbscan,length(f));
Profile.Ps_shear_k.s2 = nan(nbscan,length(f));
Profile.Ps_shear_co_k.s1 = nan(nbscan,length(f));
Profile.Ps_shear_co_k.s2 = nan(nbscan,length(f));
Profile.Pt_volt_f.t1 = nan(nbscan,length(f));
Profile.Pt_volt_f.t2 = nan(nbscan,length(f));
Profile.Pt_Tg_k.t1 = nan(nbscan,length(f));
Profile.Pt_Tg_k.t2 = nan(nbscan,length(f));
Profile.Pa_g_f.a1 = nan(nbscan,length(f));
Profile.Pa_g_f.a2 = nan(nbscan,length(f));
Profile.Pa_g_f.a3 = nan(nbscan,length(f));

Profile.Meta_Data = Meta_Data;

% ------------------------------------------------
% Loop through scans
fprintf(['Processing ' num2str(nbscan) ' scans \n'])
for p = 1:nbscan % p is the scan index.
    
    if mod(p,20)==0 && mod(p,100)~=0
        %fprintf([num2str(p) ' of ' num2str(nbscan) '\n'])
        fprintf(num2str(p))
    elseif mod(p,100)==0
        fprintf([num2str(p) '\n'])
    elseif p==nbscan
        fprintf('. \n')
    else
        fprintf('.')
    end
    
    % Get spectral data for each scan
    %     scan  =  get_scan_spectra(Profile,p);
    scan  =  get_scan_spectra_v2(Profile,p);
    
    % If there is data in the scan, add it to the profile
    if isfield(scan,'ind_ctdscan')
        Profile.ind_range_ctd(p,:) = [scan.ind_ctdscan(1),scan.ind_ctdscan(end)];
        Profile.ind_range_epsi(p,:) = [scan.ind_scan(1),scan.ind_scan(end)];
        
        %figure of merit: QC to check how the shear spectra fit panchev.
        %fom<=1 great match, fom>> 1 bad match.
        Profile.fom(p,1) = scan.fom.s1;
        Profile.fom(p,2) = scan.fom.s2;
        
        Profile.epsilon(p,1) = scan.epsilon.s1;
        Profile.epsilon(p,2) = scan.epsilon.s2;
        Profile.epsilon_co(p,1) = scan.epsilon_co.s1;
        Profile.epsilon_co(p,2) = scan.epsilon_co.s2;
        Profile.sh_fc(p,1) = scan.fc.s1;
        Profile.sh_fc(p,2) = scan.fc.s2;
        Profile.kvis(p,1) = scan.kvis;
        
        if isfield(scan.chi,'t1')
            Profile.chi(p,1) = scan.chi.t1;
            Profile.tg_fc(p,1) = scan.fc.t1;
            Profile.flag_tg_fc(p,1) = scan.flag_tg_fc.t1;
        end
        
        if isfield(scan.chi,'t2')
            Profile.chi(p,2) = scan.chi.t2;
            Profile.tg_fc(p,2) = scan.fc.t2;
            Profile.flag_tg_fc(p,2) = scan.flag_tg_fc.t2;
        end
        
        % Add spectra to profile
        Profile.k(p,:) = scan.k;
        
        Profile.Ps_volt_f.s1(p,:)     = scan.Ps_volt_f.s1(:).';
        Profile.Ps_volt_f.s2(p,:)     = scan.Ps_volt_f.s2(:).';
        Profile.Ps_shear_k.s1(p,:)    = scan.Ps_shear_k.s1(:).';
        Profile.Ps_shear_k.s2(p,:)    = scan.Ps_shear_k.s2(:).';
        Profile.Ps_shear_co_k.s1(p,:) = scan.Ps_shear_co_k.s1(:).';
        Profile.Ps_shear_co_k.s2(p,:) = scan.Ps_shear_co_k.s2(:).';
        
        Profile.Pt_volt_f.t1(p,:) = scan.Pt_volt_f.t1(:).';
        Profile.Pt_volt_f.t2(p,:) = scan.Pt_volt_f.t2(:).';
        Profile.Pt_Tg_k.t1(p,:)   = scan.Pt_Tg_k.t1(:).';
        Profile.Pt_Tg_k.t2(p,:)   = scan.Pt_Tg_k.t2(:).';
        
        Profile.Pa_g_f.a1(p,:) = scan.Pa_g_f.a1(:).';
        Profile.Pa_g_f.a2(p,:) = scan.Pa_g_f.a2(:).';
        Profile.Pa_g_f.a3(p,:) = scan.Pa_g_f.a3(:).';
        
        Profile.z(p)    = scan.z;
        Profile.w(p)    = scan.w;
        Profile.t(p)    = scan.t;
        Profile.s(p)    = scan.s;
        Profile.th(p)   = scan.th;
        Profile.sgth(p) = scan.sgth;
        
        if isfield(scan,'dnum')
            Profile.dnum(p) = scan.dnum;
        end
        
    end
    
end

Profile.Cs1a3_full = Profile.Cs1a3_full(:).';
Profile.Cs2a3_full = Profile.Cs2a3_full(:).';

%% Define varInfo and sort Profile fields
Profile = add_varInfo(Profile);
try
    Profile = sort_profile(Profile);
catch
    warning('Update sort_profile.m with the correct variable names');
end

% Save files
if saveData && isfield(Profile,'profNum')
    save_var_name = 'Profile';
    save_file_name = sprintf('Profile%03i',Profile.profNum);
    save_file = fullfile(Meta_Data.paths.profiles, ...
        [save_file_name '.mat']);
    eval(['save(''' save_file ''', ''' save_var_name ''');']);
end

% Sort Profile by standard field order
Profile = sort_profile(Profile);

end
