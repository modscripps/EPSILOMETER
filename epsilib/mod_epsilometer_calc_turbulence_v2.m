function Profile = mod_epsilometer_calc_turbulence_v2(Meta_Data,Profile_or_profNum,saveData)
% Profile = mod_epsilometer_calc_turbulence(Meta_Data,Profile_or_profNum)
%  Profile structure for Micro Structure. Inside Profile you ll find
%  temperature spectra in degC Hz^-1
%  Horizontal  velocity spectra in m^2/s^-2 Hz^-1
%  Acceleration/speed spectra in s^-1 Hz^-1
%
%  Edited summer 2020 - Nicole Couto
%  Created by Arnaud Le Boyer on 7/28/18.

if nargin<3
    saveData = 1;
end


%% Get Profile from Profile_or_profNum
if isnumeric(Profile_or_profNum) && ~isstruct(Profile_or_profNum)
    profNum = Profile_or_profNum;
    load(fullfile(Meta_Data.paths.profiles,sprintf('Profile%04.0f',profNum)));
    
elseif isstruct(Profile_or_profNum)
    Profile = Profile_or_profNum;
elseif isclassfield(Profile_or_profNum,'epsi') && isclassfield(Profile_or_profNum,'ctd') && isclassfield(Profile_or_profNum,'Meta_Data')
    Profile.Meta_Data = Profile_or_profNum.Meta_Data;
    Profile.epsi = Profile_or_profNum.epsi;
    Profile.ctd = Profile_or_profNum.ctd;
else
    error('Need epsi, ctd, and Meta_Data to calculate turbulence parameters!');
end

% Add latitude and longitude - take the mean of the first 10 seconds of the profile
% gps_sec_int = mode(seconds(days(diff(Profile_or_profNum.gps.dnum))));
%ALB I had an issue here I think this should be Profile instead of Profile_or_profNum
if isfield(Profile_or_profNum,'gps') && ~isempty(Profile_or_profNum.gps)
    gps_sec_int = mode(seconds(days(diff(Profile.gps.dnum))));
    n_sec = 10;
    idx = round(gps_sec_int*n_sec);
    nGPS = length(Profile_or_profNum.gps.latitude);
    Profile.latitude  = mean(Profile_or_profNum.gps.latitude(1:min([idx,nGPS])),'omitmissing');
    Profile.longitude = mean(Profile_or_profNum.gps.longitude(1:min([idx,nGPS])),'omitmissing');
end

%% Get epsi sampling frequency
try
    Fs_epsi = Meta_Data.AFE.FS;
catch
    Fs_epsi = Meta_Data.PROCESS.Fs_epsi;
end

%% get pitch roll rms accel

Ax=Profile.epsi.(Meta_Data.PROCESS.timeseries{1}).*0; %a1
Ay=Ax; %a3
Az=Ax; %a2

if ~isempty(Profile.epsi)
    for c=1:Meta_Data.PROCESS.nb_channels
        wh_channel=Meta_Data.PROCESS.timeseries{c};
        switch wh_channel
            case 'a1_g'
                Az=Profile.epsi.(wh_channel);
            case 'a2_g'
                Ax=Profile.epsi.(wh_channel);
            case'a3_g'
                Ay=Profile.epsi.(wh_channel);
        end
    end
end
rms_accel=rms([Ax(:);Ay(:);Az(:);]);
if sum(abs(Ax))>0 && sum(abs(Ay))>0 && sum(abs(Az))>0
    Profile.epsi.pitch = atan2 ( Ax, sqrt(Ay.^2 + Az.^2));
else
    Profile.epsi.pitch= Ax.*0;
end

if sum(abs(Ay))>0 && sum(abs(Az))>0
    Profile.epsi.roll  = atan2 ( -Ay , Az);
else
    Profile.epsi.roll= Ax.*0;
end

Profile.epsi.accel_mask=(abs(Ax)>3*rms_accel) | (abs(Ay)>3*rms_accel) |...
            (abs(Az)>3*rms_accel);



%% despiking shear and temp

%fprintf('despiking Profile\r\n')
%ALB Following the SCOR group recommendation
%I am assuming any impact between probes and stuff in the ocean is about
%50 ms = 50e-3 sec * 320Hz = 16 samples.
%I am removing any outliers that are above  3 times the standard deviation of a
%windows that is 10 * 16 samples. Arbitrarily using 10.
% I am saving the percentage of outliers in the Profile structure
movmean_window_width = Meta_Data.PROCESS.movmean_window_time*Fs_epsi;


if ~isempty(Profile.epsi)
    for c=1:Meta_Data.PROCESS.nb_channels
        wh_channel=Meta_Data.PROCESS.timeseries{c};

        Profile.epsi.([wh_channel '_raw'])=Profile.epsi.(wh_channel);
        Profile.epsi.(wh_channel)=filloutliers(Profile.epsi.(wh_channel),'linear','movmedian',movmean_window_width);

        Profile.qc.outliers.(wh_channel)= ...
            sum(Profile.epsi.(wh_channel)~= ...
            Profile.epsi.([wh_channel '_raw']))/...
            numel(Profile.epsi.(wh_channel))*100;
    end

end

%% get channels
channels = Meta_Data.PROCESS.timeseries;

nfft = Meta_Data.PROCESS.nfft;
nfftc = Meta_Data.PROCESS.nfftc;

%fpump = Meta_Data.PROCESS.ctd_fc;
%tscan = Meta_Data.PROCESS.tscan; %We don't use tscan anymore
% limit_speed = .2; % limit speed 20 cm s^{-1}
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
% NC rmoutliers only works for Matlab2018b and newer. DEV3 has Matlab 2018a
% Prmin = Meta_Data.PROCESS.Prmin(rmoutliers(Profile.ctd.P));
% Prmax = Meta_Data.PROCESS.Prmax(rmoutliers(Profile.ctd.P));
Prmin = Meta_Data.PROCESS.Prmin(Profile.ctd.P);
Prmax = Meta_Data.PROCESS.Prmax(Profile.ctd.P);

% Find ctdtime values within this pressure range
inRange = Profile.ctd.P>=Prmin & Profile.ctd.P<=Prmax;
ctdtime = Profile.ctd.time_s(inRange);

% ALB data drops can choke the code here (i.e. inRange is empty)

if sum(inRange)==0
    disp('no CTD data is that Profile')
    Profile=[];
else
% Remove epsi values outside the time period defined by ctdtime
% ALB why did I add the suffix _coh here?
% keepEpsi = Profile.epsi.time_s>=nanmin(ctdtime) & Profile.epsi.time_s<=nanmax(ctdtime);
% for iChan=1:numel(channels)
%     wh_channel = channels{iChan};
%     if ~strcmp(wh_channel,'c_count')
%         Profile_coh.(wh_channel) = Profile.epsi.(wh_channel)(keepEpsi);
%     end
% end


%% define a Pressure axis to an which I will compute epsilon and chi.
%  The spectra will be nfft long centered around P(z) +/- tscan/2.
%
Pr = ceil(min(Profile.ctd.P)):dz:floor(max(Profile.ctd.P));
nbscan = length(Pr);

%% compute coherence with a3 over the full profile.
% NC added 10/14/21 - compute coherence with a1 for earlier datasets. Check
% Meta_Data.PROCESS to determine which channel
if isfield(Meta_Data.PROCESS,'coherent_acc_chan')
    channel = [Meta_Data.PROCESS.coherent_acc_chan '_g'];
else
    channel = 'a3_g';
end


%%

% LCTD = length(Profile.ctd.P);% length of profile
% number of samples for a scan. I make sure it is always even
% N_epsi = tscan.*Fs_epsi-mod(tscan*Fs_epsi,2);
% N_ctd = tscan.*Fs_ctd-mod(tscan*Fs_ctd,2);
% Sv1 = Meta_Data.(field_name).s1.cal;
% Sv2 = Meta_Data.(field_name).s2.cal;
% Sensitivity of FPO7 probe, nominal
% dTdV(1) = Meta_Data.(field_name).t1.cal; % define in mod_epsi_temperature_spectra
% dTdV(2) = Meta_Data.(field_name).t2.cal; % define in mod_epsi_temperature_spectra
% h_freq = Meta_Data.PROCESS.h_freq;
% FPO7noise = Meta_Data.PROCESS.FPO7noise;
% % get shear probe calibration coefficient.
% if isfield(Meta_Data,'AFE')
%     field_name = 'AFE';
% elseif isfield(Meta_Data,'epsi')
%     field_name = 'epsi';
% end
% get index of the acceleration channels (useful when the number of channels is not 8)
% idxA = contains(Meta_Data.PROCESS.timeseries,'a');
% accList = Meta_Data.PROCESS.timeseries(idxA);

f = Meta_Data.PROCESS.fe;


% start creating the Profile structure
% get nb of scans in the profile
Profile.pr         =  Pr(:);
Profile.nbscan     =  nbscan;
Profile.nfft       =  nfft;
Profile.nfftc      =  nfftc;
%Profile.tscan      =  tscan; %We don't use tscan anymore
Profile.fpump      =  Meta_Data.PROCESS.ctd_fc; % arbitrary cut off frequency usually extract from coherence spectra shear/accel
Profile.f          =  f(:).';
Profile.Meta_Data  =  Meta_Data;
%% Get data for each scan and store averages in Profile

% ------------------------------------------------
% Pre-allocate profile space
Profile.dnum           = nan(nbscan,1);
Profile.z              = nan(nbscan,1);
Profile.t              = nan(nbscan,1);
Profile.w              = nan(nbscan,1);
Profile.s              = nan(nbscan,1);
Profile.th             = nan(nbscan,1);
Profile.sgth           = nan(nbscan,1);
Profile.epsilon_final  = nan(nbscan,1);
Profile.kvis           = nan(nbscan,1);
Profile.ktemp          = nan(nbscan,1);
Profile.sumPa.a1       = nan(nbscan,1);
Profile.sumPa.a2       = nan(nbscan,1);
Profile.sumPa.a3       = nan(nbscan,1);
Profile.pitch          = nan(nbscan,1);
Profile.roll           = nan(nbscan,1);
Profile.epsi_qc        = nan(nbscan,1);
Profile.qc.ratio       = nan(nbscan,2);
Profile.qc.fom         = nan(nbscan,2);
Profile.qc.a1          = nan(nbscan,1);
Profile.qc.a2          = nan(nbscan,1);

Profile.qc.a3          = nan(nbscan,1);
Profile.qc.speed       = nan(nbscan,1);
Profile.qc.pitch       = nan(nbscan,1);
Profile.qc.roll        = nan(nbscan,1);


Profile.ind_range_ctd   = nan(nbscan,2);
Profile.ind_range_epsi  = nan(nbscan,2);
Profile.fom             = nan(nbscan,2);
Profile.calib_volt      = nan(nbscan,2);
Profile.calib_vel       = nan(nbscan,2);
Profile.epsilon         = nan(nbscan,2);
Profile.epsilon_co      = nan(nbscan,2);
Profile.epsilon_method  = nan(nbscan,2);
Profile.sh_fc           = nan(nbscan,2);
Profile.sh_kc           = nan(nbscan,2);
Profile.chi             = nan(nbscan,2);
Profile.tg_fc           = nan(nbscan,2);
Profile.tg_kc           = nan(nbscan,2);
Profile.flag_tg_fc      = nan(nbscan,2);

Profile.k                = nan(nbscan,length(f));
Profile.Pt_volt_f.t1     = nan(nbscan,length(f));
Profile.Pt_volt_f.t2     = nan(nbscan,length(f));
Profile.Pt_Tg_k.t1       = nan(nbscan,length(f));
Profile.Pt_Tg_k.t2       = nan(nbscan,length(f));
Profile.Pa_g_f.a1        = nan(nbscan,length(f));
Profile.Pa_g_f.a2        = nan(nbscan,length(f));
Profile.Pa_g_f.a3        = nan(nbscan,length(f));
Profile.Ps_volt_f.s1     = nan(nbscan,length(f));
Profile.Ps_volt_f.s2     = nan(nbscan,length(f));
Profile.Ps_shear_k.s1    = nan(nbscan,length(f));
Profile.Ps_shear_k.s2    = nan(nbscan,length(f));
Profile.Ps_shear_co_k.s1 = nan(nbscan,length(f));                          %wnb shear spectrum with coherence correction
Profile.Ps_shear_co_k.s2 = nan(nbscan,length(f));                          %wnb shear spectrum with coherence correction
Profile.Cs1a3            = nan(nbscan,length(f));                          %wnb shear spectrum with coherence correction
Profile.Cs2a3            = nan(nbscan,length(f));                          %wnb shear spectrum with coherence correction

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
    
    
%==================================================================%%
%================ Core of the processing ==========================%%    
    % Get spectral data for each scan
    try
       scan  =  get_scan_spectra(Profile,p);
    catch
        disp ('pb with scan')
    end
%==================================================================%%    
%==================================================================%%    

    % If there is data in the scan, add it to the profile
    if isfield(scan,'ind_ctdscan')
        Profile.ind_range_ctd(p,:) = [scan.ind_ctdscan(1),scan.ind_ctdscan(end)];
        Profile.ind_range_epsi(p,:) = [scan.ind_scan(1),scan.ind_scan(end)];
        
        %figure of merit: QC to check how the shear spectra fit panchev.
        %fom<=1 great match, fom>> 1 bad match.
        Profile.fom(p,1) = scan.fom.s1;
        Profile.fom(p,2) = scan.fom.s2;

        Profile.calib_volt(p,1) = scan.calib_volt.s1;
        Profile.calib_volt(p,2) = scan.calib_volt.s2;
        Profile.calib_vel(p,1) = scan.calib_vel.s1;
        Profile.calib_vel(p,2) = scan.calib_vel.s2;
        
        Profile.epsilon(p,1) = scan.epsilon.s1;
        Profile.epsilon(p,2) = scan.epsilon.s2;
        Profile.epsilon_co(p,1) = scan.epsilon_co.s1;
        Profile.epsilon_co(p,2) = scan.epsilon_co.s2;
        Profile.epsilon_method(p,1) = scan.method.s1(2);
        Profile.epsilon_method(p,2) = scan.method.s2(2);
        Profile.sh_fc(p,1) = scan.fc.s1(2); %only saving the shear_co limit NC-changed to 1st index instead of 2nd. Arnaud must have changed how cutoff number is stored
        Profile.sh_fc(p,2) = scan.fc.s2(2); %only saving the shear_co limit ALB 03/06/2023 the 2nd idx is the shear_co
        Profile.sh_kc(p,1) = scan.kc.s1(2); 
        Profile.sh_kc(p,2) = scan.kc.s2(2); 
        Profile.kvis(p,1)  = scan.kvis;
        Profile.ktemp(p,1) = scan.ktemp;
        
        if isfield(scan.chi,'t1')
            Profile.chi(p,1)   = scan.chi.t1;
            Profile.tg_fc(p,1) = scan.fc.t1;
            Profile.tg_kc(p,1) = scan.kc.t1;
            Profile.flag_tg_fc(p,1) = scan.flag_tg_fc.t1;
        end
        
        if isfield(scan.chi,'t2')
            Profile.chi(p,2)   = scan.chi.t2;
            Profile.tg_fc(p,2) = scan.fc.t2;
            Profile.tg_kc(p,2) = scan.kc.t2;
            Profile.flag_tg_fc(p,2) = scan.flag_tg_fc.t2;
        end
        
        % Add spectra to profile
        Profile.k(p,:) = scan.k;
        
        Fname=fieldnames(scan.Ps_volt_f);
        for iname=1:length(Fname)
            wh_name=Fname{iname};
            Profile.Ps_volt_f.(wh_name)(p,:)     = scan.Ps_volt_f.(wh_name)(:).';
            Profile.Ps_velocity_f.(wh_name)(p,:) = scan.Ps_velocity_f.(wh_name)(:).';
            Profile.Ps_shear_k.(wh_name)(p,:)    = scan.Ps_shear_k.(wh_name)(:).';
            Profile.Ps_shear_co_k.(wh_name)(p,:) = scan.Ps_shear_co_k.(wh_name)(:).';
        end
        
        Fname=fieldnames(scan.Pt_volt_f);
        for iname=1:length(Fname)
            wh_name=Fname{iname};
            Profile.Pt_volt_f.(wh_name)(p,:) = scan.Pt_volt_f.(wh_name)(:).';
            Profile.Pt_Tg_k.(wh_name)(p,:)   = scan.Pt_Tg_k.(wh_name)(:).';
        end
        
        Fname=fieldnames(scan.Pa_g_f);
        for iname=1:length(Fname)
            wh_name=Fname{iname};
            Profile.Pa_g_f.(wh_name)(p,:) = scan.Pa_g_f.(wh_name)(:).';
            Profile.sumPa.(wh_name)(p) = scan.sumPa_g.(wh_name);
        end

        Profile.Cs1a3(p,:)=scan.Cs1a.a3(:).';
        Profile.Cs2a3(p,:)=scan.Cs2a.a3(:).';
        
        Profile.z(p)     = scan.z;
        Profile.w(p)     = scan.w;
        Profile.t(p)     = scan.t;
        Profile.s(p)     = scan.s;
        Profile.th(p)    = scan.th;
        Profile.sgth(p)  = scan.sgth;
        Profile.pitch(p) = scan.pitch;
        Profile.roll(p)  = scan.roll;

        
        if isfield(scan,'dnum')
            Profile.dnum(p) = scan.dnum;
        end
        
    end
    
end

%% add qc flag
% check on fom, sumPa, calib_vel ...

% sumPa
rmsa1=rms(log10(Profile.sumPa.a1),'omitnan');
stda1=std(log10(Profile.sumPa.a1),'omitnan');
rmsa2=rms(log10(Profile.sumPa.a2),'omitnan');
stda2=std(log10(Profile.sumPa.a2),'omitnan');
rmsa3=rms(log10(Profile.sumPa.a3),'omitnan');
stda3=std(log10(Profile.sumPa.a3),'omitnan');

qc.a1=log10(Profile.sumPa.a1)>(-rmsa1+2.*stda1);
qc.a2=log10(Profile.sumPa.a2)>(-rmsa2+2.*stda2);
qc.a3=log10(Profile.sumPa.a3)>(-rmsa3+2.*stda3);

% figure of merit
qc.fom=Profile.fom>2;

%pitch
qc.pitch = (Profile.pitch-mean(Profile.pitch,'omitnan'))>5;
qc.roll  = (Profile.roll-mean(Profile.roll,'omitnan'))>5;

% speed
try
    qc.speed = Profile.w>Meta_Data.PROCESS.speed_limit;
catch
    disp('No speed limit in Meta_Data.PROCESS')
    Meta_Data.PROCESS.speed_limit=.2;
    qc.speed = Profile.w<Meta_Data.PROCESS.speed_limit;
end

%epsi1/epsi2 ratio
qc.ratio(:,1) = (Profile.epsilon_co(:,1)./Profile.epsilon_co(:,2))>3;
qc.ratio(:,2) = (Profile.epsilon_co(:,2)./Profile.epsilon_co(:,1))>3;

final_epsi=Profile.epsilon_co;
final_epsi(qc.fom)=nan;

Profile.epsilon_final=mean(final_epsi,2,'omitnan');
Profile.epsilon_final(qc.ratio(:,1))=final_epsi(qc.ratio(:,1),2);
Profile.epsilon_final(qc.ratio(:,2))=final_epsi(qc.ratio(:,2),1);



qc_array = [qc.a3 qc.a2 qc.a1 qc.roll qc.pitch qc.speed];
qc_bit   =  bin2dec(num2str(qc_array));
Profile.epsi_qc = qc_bit;

Profile.epsilon_final(qc_bit>0) =nan;

qc_array = [qc.ratio qc.fom qc.a3 qc.a2 qc.a1 qc.roll qc.pitch qc.speed];
qc_bit   =  bin2dec(num2str(qc_array));
Profile.epsi_qc = qc_bit;
Profile.qc      = qc;



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
    save_file_name = sprintf('Profile%04i',Profile.profNum);
    save_file = fullfile(Meta_Data.paths.profiles, ...
        [save_file_name '.mat']);
    eval(['save(''' save_file ''', ''' save_var_name ''');']);
end

% Sort Profile by standard field order
Profile = sort_profile(Profile);

end
end
