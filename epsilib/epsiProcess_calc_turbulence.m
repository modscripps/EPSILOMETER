function data_struct = epsiProcess_calc_turbulence(Meta_Data,data_struct,saveData)
% Profile = mod_epsilometer_calc_turbulence(Meta_Data,Profile_or_profNum)
%  Profile structure for Micro Structure. Inside Profile you ll find
%  temperature spectra in degC Hz^-1
%  Horizontal  velocity spectra in m^2/s^-2 Hz^-1
%  Acceleration/speed spectra in s^-1 Hz^-1
%
%  Created by Arnaud Le Boyer on 7/28/18.
%  Edited summer 2020 - Nicole Couto
% NC - May 2021 - copied from mod_epsilometer_calc_turbulence_v2.m and
% edited to calculate turbulence variables at intervals of time rather than
% pressure

if nargin<3
    saveData = 1;
end

%% Get epsi sampling frequency
try
    Fs_epsi = Meta_Data.AFE.FS;
catch
    Fs_epsi = Meta_Data.PROCESS.Fs_epsi;
end

%% despiking shear and temp

%fprintf('despiking Profile\r\n')
%ALB Following the SCOR group recommendation
%I am assuming any impact between probes and stuff in the ocean is about
%50 ms = 50e-3 sec * 320Hz = 16 samples.
%I am removing any outliers that are above  3 times the standard deviation of a
%windows that is 10 * 16 samples. Arbitrarily using 10.
% I am saving the percentage of outliers in the Profile structure
movmean_window_width = Meta_Data.PROCESS.movmean_window_time*Fs_epsi;
%movmean_window_width=10*16;
if ~isempty(data_struct.epsi)
    for c=1:Meta_Data.PROCESS.nb_channels
        wh_channel=Meta_Data.PROCESS.timeseries{c};
        if ~strcmp(wh_channel,'c_count')
            data_struct.epsi.([wh_channel '_raw'])=data_struct.epsi.(wh_channel);
            data_struct.epsi.(wh_channel)=  ...
                filloutliers(data_struct.epsi.([wh_channel '_raw']), ...
                'linear','movmean',movmean_window_width);
            data_struct.qc.outliers.(wh_channel)= ...
                sum(data_struct.epsi.(wh_channel)~= ...
                data_struct.epsi.([wh_channel '_raw']))/...
                numel(data_struct.epsi.(wh_channel))*100;
        end
    end
    
end

%% get channels
channels = Meta_Data.PROCESS.timeseries;

nfft = Meta_Data.PROCESS.nfft;
nfftc = Meta_Data.PROCESS.nfftc;

%fpump = Meta_Data.PROCESS.ctd_fc;
%tscan = Meta_Data.PROCESS.tscan; %We don't use tscan anymore
% limit_speed = .2; % limit speed 20 cm s^{-1}
dz  = Meta_Data.PROCESS.dz;
dt  = Meta_Data.PROCESS.dt;

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
% % Find max and min pressure
% Prmin = Meta_Data.PROCESS.Prmin(data_struct.ctd.P);
% Prmax = Meta_Data.PROCESS.Prmax(data_struct.ctd.P);
% 
% % Find ctdtime values within this pressure range
% inRange = data_struct.ctd.P>=Prmin & data_struct.ctd.P<=Prmax;
% ctdtime = data_struct.ctd.time_s(inRange);
% 
% 
% % Remove epsi values outside the time period defined by ctdtime
% %ALB why did I add the suffix _coh here?
% keepEpsi = data_struct.epsi.time_s>=nanmin(ctdtime) & data_struct.epsi.time_s<=nanmax(ctdtime);
% for iChan=1:numel(channels)
%     wh_channel = channels{iChan};
%     if ~strcmp(wh_channel,'c_count')
%         Profile_coh.(wh_channel) = data_struct.epsi.(wh_channel)(keepEpsi);
%     end
% end
% 
% if isfield(Meta_Data.PROCESS, "multivariate")
%     if (Meta_Data.PROCESS.multivariate)
%         for iChan=1:numel(channels)
%             wh_channel = channels{iChan};
%             if ~strcmp(wh_channel,'c_count')
%                 Profile_coh.(wh_channel) = data_struct.epsi.(wh_channel);
%             end
%         end
%         %ALB add ctd.P and ctd.time_s so I can compute multivariate coherence on
%         %the same grid as the futur scans
%         Profile_coh.time_s=data_struct.epsi.time_s;
%         Profile_coh.ctd.P=data_struct.ctd.P;
%         Profile_coh.ctd.time_s=data_struct.ctd.time_s;
%     end
% end

%% define a time axis to an which I will compute epsilon and chi.
%  The spectra will be nfft long centered around t(idx) +/- tscan/2.
%
t_array = ceil(min(data_struct.ctd.time_s)):dt:floor(max(data_struct.ctd.time_s));
nbscan = length(t_array);

 %% compute coherence with a3 over the full data_struct.
% % NC added 10/14/21 - compute coherence with a1 for earlier datasets. Check
% % Meta_Data.PROCESS to determine which channel
% if isfield(Meta_Data.PROCESS,'coherent_acc_chan')
%     channel = [Meta_Data.PROCESS.coherent_acc_chan '_g'];
% else
%     channel = 'a3_g';
% end
% 
% if isfield(Meta_Data.PROCESS, "multivariate")
%     % if the user asked to compute the multivariate and if it is not
%     % already done
%     if (Meta_Data.PROCESS.multivariate) && ~isfield(Profile,"Cs1a_mv") 
%         
%         [data_struct.Cs1a_mv,data_struct.Cs2a_mv] = ...
%             mod_efe_profile_multivariate_coherence(Profile_coh,Pr,Meta_Data);
%         % save profile cause it takes for eve to compute mv.
%         save_var_name = 'Profile';
%         save_file_name = sprintf('Profile%03i',data_struct.profNum);
%         save_file = fullfile(Meta_Data.paths.profiles, ...
%             [save_file_name '.mat']);
%         eval(['save(''' save_file ''', ''' save_var_name ''');']);
%     end
% end
% 
% [data_struct.Cs1a3_full,data_struct.Cs2a3_full,...
%     ~,~,~] = mod_efe_scan_coherence(Profile_coh,channel,Meta_Data);

f = Meta_Data.PROCESS.fe;


% start creating the Profile structure
% get nb of scans in the profile
data_struct.time_s     =  t_array(:);
data_struct.nbscan     =  nbscan;
data_struct.nfft       =  nfft;
data_struct.nfftc      =  nfftc;
data_struct.fpump      =  Meta_Data.PROCESS.ctd_fc; % arbitrary cut off frequency usually extract from coherence spectra shear/accel
data_struct.f          =  f(:).';

%% Get data for each scan and store averages in Profile

% ------------------------------------------------
% Pre-allocate profile space
data_struct.dnum           = nan(nbscan,1);
data_struct.z              = nan(nbscan,1);
data_struct.t              = nan(nbscan,1);
data_struct.w              = nan(nbscan,1);
data_struct.s              = nan(nbscan,1);
data_struct.th             = nan(nbscan,1);
data_struct.sgth           = nan(nbscan,1);
data_struct.epsilon_final  = nan(nbscan,1);
data_struct.kvis           = nan(nbscan,1);
data_struct.epsi_qc_final  = nan(nbscan,1);

data_struct.ind_range_ctd   = nan(nbscan,2);
data_struct.ind_range_epsi  = nan(nbscan,2);
data_struct.fom             = nan(nbscan,2);
data_struct.calib_volt      = nan(nbscan,2);
data_struct.calib_vel       = nan(nbscan,2);
data_struct.epsilon         = nan(nbscan,2);
data_struct.epsilon_co      = nan(nbscan,2);
data_struct.epsilon_mv      = nan(nbscan,2);
data_struct.sh_fc           = nan(nbscan,2);
data_struct.chi             = nan(nbscan,2);
data_struct.tg_fc           = nan(nbscan,2);
data_struct.flag_tg_fc      = nan(nbscan,2);
data_struct.epsi_qc         = nan(nbscan,2);

data_struct.k                = nan(nbscan,length(f));
data_struct.Pt_volt_f.t1     = nan(nbscan,length(f));
data_struct.Pt_volt_f.t2     = nan(nbscan,length(f));
data_struct.Pt_Tg_k.t1       = nan(nbscan,length(f));
data_struct.Pt_Tg_k.t2       = nan(nbscan,length(f));
data_struct.Pa_g_f.a1        = nan(nbscan,length(f));
data_struct.Pa_g_f.a2        = nan(nbscan,length(f));
data_struct.Pa_g_f.a3        = nan(nbscan,length(f));
data_struct.Ps_volt_f.s1     = nan(nbscan,length(f));
data_struct.Ps_volt_f.s2     = nan(nbscan,length(f));
data_struct.Ps_shear_k.s1    = nan(nbscan,length(f));
data_struct.Ps_shear_k.s2    = nan(nbscan,length(f));
data_struct.Ps_shear_co_k.s1 = nan(nbscan,length(f));                          %wnb shear spectrum with coherence correction
data_struct.Ps_shear_co_k.s2 = nan(nbscan,length(f));                          %wnb shear spectrum with coherence correction
data_struct.Ps_shear_mv_k.s1 = nan(nbscan,length(f));                          %wnb shear spectrum with multivariate correction
data_struct.Ps_shear_mv_k.s2 = nan(nbscan,length(f));                          %wnb shear spectrum with multivariate correction

data_struct.Meta_Data = Meta_Data;

% ------------------------------------------------
% Loop through scans
fprintf(['Processing ' num2str(nbscan) ' scans \n'])
for s = 1:nbscan % s is the scan index.
    
    if mod(s,20)==0 && mod(s,100)~=0
        %fprintf([num2str(p) ' of ' num2str(nbscan) '\n'])
        fprintf(num2str(s))
    elseif mod(s,100)==0
        fprintf([num2str(s) '\n'])
    elseif s==nbscan
        fprintf('. \n')
    else
        fprintf('.')
    end
    
    
%==================================================================%%
%================ Core of the processing ==========================%%    
    % Get spectral data for each scan
    scan  =  get_scan_spectra_t(data_struct,s);
%==================================================================%%    
%==================================================================%%    

    % If there is data in the scan, add it to the profile
    if isfield(scan,'ind_ctdscan')
        data_struct.ind_range_ctd(s,:) = [scan.ind_ctdscan(1),scan.ind_ctdscan(end)];
        data_struct.ind_range_epsi(s,:) = [scan.ind_scan(1),scan.ind_scan(end)];
        
        %figure of merit: QC to check how the shear spectra fit panchev.
        %fom<=1 great match, fom>> 1 bad match.
        data_struct.fom(s,1) = scan.fom.s1;
        data_struct.fom(s,2) = scan.fom.s2;

        data_struct.calib_volt(s,1) = scan.calib_volt.s1;
        data_struct.calib_volt(s,2) = scan.calib_volt.s2;
        data_struct.calib_vel(s,1) = scan.calib_vel.s1;
        data_struct.calib_vel(s,2) = scan.calib_vel.s2;
        
        data_struct.epsilon(s,1) = scan.epsilon.s1;
        data_struct.epsilon(s,2) = scan.epsilon.s2;
        data_struct.epsilon_co(s,1) = scan.epsilon_co.s1;
        data_struct.epsilon_co(s,2) = scan.epsilon_co.s2;
        data_struct.epsilon_final(s) = scan.epsilon_final;
        data_struct.epsi_qc(s,1) = scan.epsi_qc(1);
        data_struct.epsi_qc(s,2) = scan.epsi_qc(2);
        data_struct.epsi_qc_final(s) = scan.epsi_qc_final;
        data_struct.sh_fc(s,1) = scan.fc.s1(2); %only saving the shear_co limit
        data_struct.sh_fc(s,2) = scan.fc.s2(2); %only saving the shear_co limit
        data_struct.kvis(s,1) = scan.kvis;
        
        if isfield(scan.chi,'t1')
            data_struct.chi(s,1) = scan.chi.t1;
            data_struct.tg_fc(s,1) = scan.fc.t1;
            data_struct.flag_tg_fc(s,1) = scan.flag_tg_fc.t1;
        end
        
        if isfield(scan.chi,'t2')
            data_struct.chi(s,2) = scan.chi.t2;
            data_struct.tg_fc(s,2) = scan.fc.t2;
            data_struct.flag_tg_fc(s,2) = scan.flag_tg_fc.t2;
        end
        
        % Add spectra to profile
        data_struct.k(s,:) = scan.k;
        
        data_struct.Ps_volt_f.s1(s,:)     = scan.Ps_volt_f.s1(:).';
        data_struct.Ps_volt_f.s2(s,:)     = scan.Ps_volt_f.s2(:).';
        data_struct.Ps_velocity_f.s1(s,:)     = scan.Ps_velocity_f.s1(:).';
        data_struct.Ps_velocity_f.s2(s,:)     = scan.Ps_velocity_f.s2(:).';
        data_struct.Ps_shear_k.s1(s,:)    = scan.Ps_shear_k.s1(:).';
        data_struct.Ps_shear_k.s2(s,:)    = scan.Ps_shear_k.s2(:).';
        data_struct.Ps_shear_co_k.s1(s,:) = scan.Ps_shear_co_k.s1(:).';
        data_struct.Ps_shear_co_k.s2(s,:) = scan.Ps_shear_co_k.s2(:).';
        
        data_struct.Pt_volt_f.t1(s,:) = scan.Pt_volt_f.t1(:).';
        data_struct.Pt_volt_f.t2(s,:) = scan.Pt_volt_f.t2(:).';
        data_struct.Pt_Tg_k.t1(s,:)   = scan.Pt_Tg_k.t1(:).';
        data_struct.Pt_Tg_k.t2(s,:)   = scan.Pt_Tg_k.t2(:).';
        
        data_struct.Pa_g_f.a1(s,:) = scan.Pa_g_f.a1(:).';
        data_struct.Pa_g_f.a2(s,:) = scan.Pa_g_f.a2(:).';
        data_struct.Pa_g_f.a3(s,:) = scan.Pa_g_f.a3(:).';
        
        data_struct.z(s)    = scan.z;
        data_struct.w(s)    = scan.w;
        data_struct.t(s)    = scan.t;
        data_struct.s(s)    = scan.s;
        data_struct.th(s)   = scan.th;
        data_struct.sgth(s) = scan.sgth;
        
        if isfield(scan,'dnum')
            data_struct.dnum(s) = scan.dnum;
        end
        
    end
    
end

if isfield(data_struct,'Cs1a3_full')
    data_struct.Cs1a3_full = data_struct.Cs1a3_full(:).';
    data_struct.Cs2a3_full = data_struct.Cs2a3_full(:).';
end

%% Define varInfo and sort Profile fields
data_struct = add_varInfo(data_struct);
try
    data_struct = sort_profile(data_struct);
catch
    warning('Update sort_profile.m with the correct variable names');
end

% % Save files
% if saveData && isfield(Profile,'profNum')
%     save_var_name = 'Profile';
%     save_file_name = sprintf('Profile%03i',data_struct.profNum);
%     save_file = fullfile(Meta_Data.paths.profiles, ...
%         [save_file_name '.mat']);
%     eval(['save(''' save_file ''', ''' save_var_name ''');']);
% end

% Sort Profile by standard field order
%data_struct = sort_profile(data_struct);

end
