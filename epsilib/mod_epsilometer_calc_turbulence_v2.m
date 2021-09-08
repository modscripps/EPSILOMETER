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
    load(fullfile(Meta_Data.L1path,sprintf('Profile%03.0f',profNum)));
    %eval(['Profile = ' sprintf('Profile%03.0f',profNum) ';']);
elseif isstruct(Profile_or_profNum)
    Profile = Profile_or_profNum;
elseif isclassfield(Profile_or_profNum,'epsi') && isclassfield(Profile_or_profNum,'ctd') && isclassfield(Profile_or_profNum,'Meta_Data')
    Profile.Meta_Data = Profile_or_profNum.Meta_Data;
    Profile.epsi = Profile_or_profNum.epsi;
    Profile.ctd = Profile_or_profNum.ctd;
else
    error('Need epsi, ctd, and Meta_Data to calculate turbulence parameters!');
end


%% get channels
channels = Meta_Data.PROCESS.timeseries;

nfft = Meta_Data.PROCESS.nfft;
nfftc = Meta_Data.PROCESS.nfftc;
Fs_epsi = Meta_Data.AFE.FS;
Fs_ctd = Meta_Data.PROCESS.Fs_ctd;

fpump = Meta_Data.PROCESS.ctd_fc;
tscan = Meta_Data.PROCESS.tscan;
limit_speed = .2; % limit speed 20 cm s^{-1}
dz  =   Meta_Data.PROCESS.dz;

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
    Profile_coh.(wh_channel) = Profile.epsi.(wh_channel)(keepEpsi);
end

% compute coherence with a3 over the full profile.
[Profile.Cs1a3_full,Profile.Cs2a3_full,...
        ~,~,~] = mod_efe_scan_coherence(Profile_coh,'a3_g',Meta_Data);

% NC - Profile already has dPdt
% %% Get dPdt
% switch Meta_Data.vehicle_name
%     case 'FISH'
%         Profile.dPdt  =  compute_fallrate_downcast(Profile);
%     case {'WW','Seacycler'}
%         % TODO: make the P from the WW CTD in the same unit as SEABIRD
%         Profile.dPdt  =  compute_speed_upcast(Profile);
%         Profile.dPdt  =  -Profile.dPdt/1e7;
% end

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
Sv1 = Meta_Data.AFE.s1.cal;
Sv2 = Meta_Data.AFE.s2.cal;
% Sensitivity of FPO7 probe, nominal
dTdV(1) = Meta_Data.AFE.t1.cal; % define in mod_epsi_temperature_spectra
dTdV(2) = Meta_Data.AFE.t2.cal; % define in mod_epsi_temperature_spectra

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
Profile.t = nan(nbscan,1);
Profile.w = nan(nbscan,1);
Profile.s = nan(nbscan,1);
Profile.dnum = nan(nbscan,1);
Profile.ind_range_ctd = nan(nbscan,2);
Profile.ind_range_epsi = nan(nbscan,2);
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
    scan  =  get_scan_spectra(Profile,p);
    
    % If there is data in the scan, add it to the profile
    if isfield(scan,'ind_ctdscan')
        Profile.ind_range_ctd(p,:) = [scan.ind_ctdscan(1),scan.ind_ctdscan(end)];
        Profile.ind_range_epsi(p,:) = [scan.ind_scan(1),scan.ind_scan(end)];
        
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
        
        Profile.Ps_volt_f.s1(p,:) = scan.Ps_volt_f.s1(:).';
        Profile.Ps_volt_f.s2(p,:) = scan.Ps_volt_f.s2(:).';
        Profile.Ps_shear_k.s1(p,:) = scan.Ps_shear_k.s1(:).';
        Profile.Ps_shear_k.s2(p,:) = scan.Ps_shear_k.s2(:).';
        Profile.Ps_shear_co_k.s1(p,:) = scan.Ps_shear_co_k.s1(:).';
        Profile.Ps_shear_co_k.s2(p,:) = scan.Ps_shear_co_k.s2(:).';       
        
        Profile.Pt_volt_f.t1(p,:) = scan.Pt_volt_f.t1(:).';
        Profile.Pt_volt_f.t2(p,:) = scan.Pt_volt_f.t2(:).';
        Profile.Pt_Tg_k.t1(p,:) = scan.Pt_Tg_k.t1(:).';
        Profile.Pt_Tg_k.t2(p,:) = scan.Pt_Tg_k.t2(:).';
        
        Profile.Pa_g_f.a1(p,:) = scan.Pa_g_f.a1(:).';
        Profile.Pa_g_f.a2(p,:) = scan.Pa_g_f.a2(:).';
        Profile.Pa_g_f.a3(p,:) = scan.Pa_g_f.a3(:).';

        Profile.w(p) = scan.w;
        Profile.t(p) = scan.t;
        Profile.s(p) = scan.s;
        if isfield(scan,'dnum')
            Profile.dnum(p) = scan.dnum;
        end
        
        % NC 8/25/21 - commented out because these fields already exist
        % inside epsi field
%         % Fill in the epsi channels
%         for c = 1:length(channels)
%             wh_channel = channels{c};   
%             Profile.(wh_channel)(scan.ind_scan) = scan.(wh_channel);
%         end
        
     
    end
    
end

%% Define varInfo and sort Profile fields
% The order of the fields in varInfo will define the order of fields in
% Profile, so make sure there are the same number of fields!
if isfield(Profile.ctd,'dnum')
    Profile.varInfo.dnum = {'datenum','Matlab datenum'};
    Profile.varInfo.time_s = {'time','seconds since Jan 1 1970'};
else
    Profile.varInfo.time_s = {'time','seconds since power on'};
end
Profile.varInfo.P = {'CTD P','db'};
Profile.varInfo.dPdt = {'CTD diff(P)/diff(ctdtime)','db s^{-1}'};
Profile.varInfo.T = {'CTD temperature','C'};
Profile.varInfo.C = {'CTD conductivity',''};
Profile.varInfo.S = {'CTD salinity','psu'};
Profile.varInfo.sig = {'CTD potential density (sigma-theta)',''};
%Profile.varInfo.EPSInbsample = {'',''};
%Profile.varInfo.epsitime = {'Epsi time','Matlab datenum'};
Profile.varInfo.a1_g = {'acceleration sensor 1 timeseries in this scan','[g]'};
Profile.varInfo.a2_g = {'acceleration sensor 2 timeseries in this scan','[g]'};
Profile.varInfo.a3_g = {'acceleration sensor 3 timeseries in this scan','[g]'};
Profile.varInfo.s1_volt = {'shear sensor 1 timeseries in this scan','Volts'};
Profile.varInfo.s2_volt = {'shear sensor 2 timeseries in this scan','Volts'};
Profile.varInfo.t1_volt = {'temperature sensor 1 timeseries in this scan','Volts'};
Profile.varInfo.t2_volt = {'temperature sensor 2 timeseries in this scan','Volts'};
if isfield(Meta_Data.PROCESS.timeseries,'c_count')
Profile.varInfo.c_count = {'additional sensor timseries in the scan' 'count'};
end
Profile.varInfo.nbscan = {'',''};
Profile.varInfo.nfft = {'',''};
Profile.varInfo.nfftc = {'',''};
Profile.varInfo.tscan = {'length of scan window','s'};
Profile.varInfo.fpump = {'',''};
Profile.varInfo.dnum = {'datenum','Matlab datenum'}; 
Profile.varInfo.pr = {'CTD pressure','db'};  
Profile.varInfo.w = {'fall speed','db s^{-1}'};
Profile.varInfo.t = {'temperature','C'};
Profile.varInfo.s = {'salinity','psu'};  
Profile.varInfo.kvis = {'kinematic viscosity',''};  
Profile.varInfo.epsilon = {'turbulent kinetic energy dissipation rate calculated from Ps_shear_k', ''};
Profile.varInfo.epsilon_co = {'turbulent kinetic energy dissipation rate calculated from Ps_shear_co_k', ''};
Profile.varInfo.chi = {'temperature gradient dissipation rate',''};
Profile.varInfo.sh_fc = {'shear cutoff frequency, 1=uncorrected, 2=coherence-corrected', 'Hz'};
Profile.varInfo.tg_fc = {'temperature gradient cutoff frequency, 1=uncorrected, 2=coherence-corrected', 'Hz'};
Profile.varInfo.flag_tg_fc = {'temperature gradient cut off frequency is very high','0/1'};
Profile.varInfo.ind_range_ctd = {'1st and last indices of CTD arrays in this Profile that go into each scan window',''};
Profile.varInfo.ind_range_epsi = {'1st and last indices of Epsi arrays in this Profile that go into each scan window',''};
Profile.varInfo.f = {'frequency','Hz'};
Profile.varInfo.k = {'wavenumber','cycles m^-^1'};
Profile.varInfo.Pa_g_f = {'accleration frequency power spectrum', '[g]^2 Hz^{-1}'};
Profile.varInfo.Ps_volt_f = {'shear frequency power spectrum', 'Volts^2 Hz^{-1}'};
Profile.varInfo.Ps_shear_k = {'shear wavenumber power spectrum', 's{-1} cpm^{-1}'};
Profile.varInfo.Ps_shear_co_k = {'coherence-corrected shear frequency power spectrum (full profile coherence with a3 channel has been removed)', ''};
Profile.varInfo.Pt_volt_f = {'temperature frequency power spectrum','Volts^2 Hz{-1}'};
Profile.varInfo.Pt_Tg_k = {'temperature gradient wavenumber power spectrum', 'C^2 s{-1} cpm^{-1}'};
Profile.varInfo.Cs1a3_full = {'coherence betwen s1 and a3 channels between Meta_Data.PROCESS.Prmin and Meta_Data.PROCESS.Prmax',''};
Profile.varInfo.Cs2a3_full = {'coherence betwen s2 and a3 channels between Meta_Data.PROCESS.Prmin and Meta_Data.PROCESS.Prmax',''};

% Profile fields to remove
fieldsToRemove = {'P_raw','T_raw','S_raw','C_raw','PT_raw','timestamp',...
't1_count','t2_count','s1_count','s2_count','a1_count','a2_count','a3_count'};
for iField=1:length(fieldsToRemove)
   if isfield(Profile,fieldsToRemove{iField})
       Profile = rmfield(Profile,fieldsToRemove{iField});
   end
end


% Sort Profile field names by order of varInfo as written above
varFields = fields(Profile.varInfo);
profFields = fields(Profile);

% Save files
if saveData && isfield(Profile,'profNum')
    save_var_name = 'Profile';
    save_file_name = sprintf('Profile%03i',Profile.profNum);
    save_file = fullfile(Meta_Data.L1path, ...
                        [save_file_name '.mat']);
    eval(['save(''' save_file ''', ''' save_var_name ''');']);
end


% 

% try
%     Profile = orderfields(Profile,['Meta_Data';'varInfo';varFields]);
% catch
%     error('There are more Profile fields than are defined in varInfo')
% end
% 




end
