function [Profile] = add_avg_psd_to_profile(Profile)
% Add average values of power spectral densities and coherences between 10-45 Hz
%
% Nicole Couto | October 2020
% -------------------------------------------------------------------------

% Get frequency arrays and define max and min f
try
    f = Profile.f;
catch
    disp('no Profile.f')
end
f1 = 10;
f2 = 45;

% Get constants for computing coherences
nfft = Profile.Meta_Data.PROCESS.nfft;
nfftc = Profile.Meta_Data.PROCESS.nfftc;
try
    Fs = Profile.Meta_Data.PROCESS.Fs_epsi;
catch
    Fs = Profile.Meta_Data.AFE.FS;
end

% Average power spectral densities
% -------------------------------------------------------------------------
% Get spectrum fields
psdList = {'Ps_volt_f','Pt_volt_f','Pa_g_f'};
for iP=1:numel(psdList)
    
    % List channels within this field
    channels = fields(Profile.(psdList{iP}));
    
    for iC=1:numel(channels)
        
        % Get current channel
        wh_channel = channels{iC};
        
        % Loop through all scan windows
        for iS=1:Profile.nbscan
            
            % Get current power spectrum
            psd = Profile.(psdList{iP}).(wh_channel)(iS,:);
            
            % Find average power spectral density between f1 and f2
            psdAvg = avg_psd_in_frange(psd,f,f1,f2);
            
            % Save average psd to Profile
            fieldStr = [psdList{iP}(1:end-2),'_10_45Hz_avg'];
            Profile.(fieldStr).(wh_channel)(iS,1) = psdAvg;
            
        end
    end
end

% Average power spectral densities
% -------------------------------------------------------------------------
% Get coherences between all shear/temperature/acceleration sensors and all acceleration
% sensors
combos = {'a1','s1','Pa_g_f','Ps_volt_f';...
    'a1','s2','Pa_g_f','Ps_volt_f';...
    'a1','t1','Pa_g_f','Pt_volt_f';...
    'a1','t2','Pa_g_f','Pt_volt_f';...
    'a1','a2','Pa_g_f','Pa_g_f';...
    'a1','a3','Pa_g_f','Pa_g_f';...
    'a2','s1','Pa_g_f','Ps_volt_f';...
    'a2','s2','Pa_g_f','Ps_volt_f';...
    'a2','t1','Pa_g_f','Pt_volt_f';...
    'a2','t2','Pa_g_f','Pt_volt_f';...
    'a2','a3','Pa_g_f','Pa_g_f';...
    'a3','s1','Pa_g_f','Ps_volt_f';...
    'a3','s2','Pa_g_f','Ps_volt_f';...
    'a3','t1','Pa_g_f','Pt_volt_f';...
    'a3','t2','Pa_g_f','Pt_volt_f'};

for iC=1:length(combos)
    
    % Get current combination name and two power spectra
    comboName = [combos{iC,1},combos{iC,2}];
    psd1 = Profile.(combos{iC,3}).(combos{iC,1});
    psd2 = Profile.(combos{iC,4}).(combos{iC,2});
    
    % Loop through all scan windows
    for iS=1:Profile.nbscan
        
        % Compute the coherence if there are at least more than 8 samples
        if sum(~isnan(detrend(psd1(iS,:))))>8 && sum(~isnan(detrend(psd2(iS,:))))>8
            Coh = mscohere(detrend(psd1(iS,:)),detrend(psd2(iS,:)),nfftc,[],nfft,Fs);
            
            % Find average coherence
            Cavg = avg_psd_in_frange(Coh,f,f1,f2);
            
            % Save to Profile with comboName
            Profile.Coh_10_45Hz_avg.(comboName)(iS,1) = Cavg;
            
        else
            
            % If there aren't enough samples to compute coherence, make it
            % nan
            Profile.Coh_10_45Hz_avg.(comboName)(iS,1) = nan;
            
        end  
    end   
end