function Meta_Data=load_metadata_from_deployment_log(path_mission)

path_mission = fileparts(logFile);

% Meta_Data.process=process_dir;
% Meta_Data.path_mission=path_mission;
% 
% fid=fopen(filename);
% header=fgetl(fid);
% 
% % read filename
% disp(['Read ' filename])
% nb_line=1;
% while(feof(fid)==0)
%     A{nb_line}=fgetl(fid);
%     nb_line=nb_line+1;
% end
% fclose(fid);
% disp(['Done reading ' filename])
% % done reading filename
% 
% % define Meta_Data parameters
% mission=strsplit(A{1},',');
% Meta_Data.mission=mission{2};
% 
% vehicle=strsplit(A{3},',');
% Meta_Data.vehicle=vehicle{2};
% 
% vehicle_name=strsplit(A{4},',');
% Meta_Data.vehicle_name=vehicle_name{2};
% 
% 
% deployment=strsplit(A{5},',');
% Meta_Data.deployment=deployment{2};
% 
% NFFT=strsplit(A{6},',');
% Meta_Data.PROCESS.nfft=str2double(NFFT{2});
% 
% Fs_epsi=strsplit(A{7},',');
% Meta_Data.PROCESS.Fs_epsi=str2double(Fs_epsi{2});
% 
% s1_SN=strsplit(A{8},',');
% Meta_Data.epsi.s1.SN=s1_SN{2};
% 
% s2_SN=strsplit(A{9},',');
% Meta_Data.epsi.s2.SN=s2_SN{2};
% 
% 
% t1_SN=strsplit(A{10},',');
% Meta_Data.epsi.t1.SN=t1_SN{2};
% 
% t2_SN=strsplit(A{11},',');
% Meta_Data.epsi.t2.SN=t2_SN{2};
% 
% CTD_name=strsplit(A{12},',');
% Meta_Data.aux1.name = CTD_name{2};
% 
% CTD_SN=strsplit(A{13},',');
% Meta_Data.aux1.SN   = sprintf('%04.0f',str2num(CTD_SN{2})); %Make 4-digits just in case it's not in the log file
% 
% channels=strsplit(A{14},',');
% channels=channels(2:end-1);
% channels{1}=channels{1}(2:end);
% channels{end}=channels{end}(1:end-1);
% Meta_Data.PROCESS.channels=channels;
% % number of channels
% Meta_Data.PROCESS.nb_channels=length(Meta_Data.PROCESS.channels);
% 
% % recording mode (STREAMING or SD)
% recording_mode=strsplit(A{15},',');
% Meta_Data.PROCESS.recording_mod=recording_mode{2};
% 
% starttime=strsplit(A{16},',');
% Meta_Data.starttime=datenum(starttime{2});
% 
% 
% MADRE_rev=strsplit(A{17},',');
% Meta_Data.MADRE.rev=MADRE_rev{2};
% MADRE_SN=strsplit(A{18},',');
% Meta_Data.MADRE.SN=MADRE_SN{2};
% 
% MAP_rev=strsplit(A{19},',');
% Meta_Data.MAP.rev=MAP_rev{2};
% 
% MAP_SN=strsplit(A{20},',');
% Meta_Data.MAP.SN=MAP_SN{2};
% 
% MAP_temp=strsplit(A{21},',');
% Meta_Data.MAP.temperature=MAP_temp{2};
% MAP_shear=strsplit(A{22},',');
% Meta_Data.MAP.shear=MAP_shear{2};
% 
% 
% Firm_version=strsplit(A{23},',');
% Meta_Data.Firmware.version=Firm_version{2};
% Firm_ADCshear=strsplit(A{24},',');
% Meta_Data.Firmware.ADCshear=Firm_ADCshear{2};
% Firm_ADCFPO7=strsplit(A{25},',');
% Meta_Data.Firmware.ADC_FPO7=Firm_ADCFPO7{2}';
% Firm_ADCaccell=strsplit(A{26},',');
% Meta_Data.Firmware.ADC_accellerometer=Firm_ADCaccell{2};
% Firm_ADCcond=strsplit(A{27},',');
% Meta_Data.Firmware.ADC_cond=Firm_ADCcond{2};
% Firm_ADCfilter=strsplit(A{28},',');
% Meta_Data.Firmware.ADCfilter=Firm_ADCfilter{2};
% % done define Meta_Data parameters
% 
% % define the frequency axes for the spectral computation
% [~,Meta_Data.PROCESS.fe] = pwelch(0*(1:Meta_Data.PROCESS.nfft),...
%                 Meta_Data.PROCESS.nfft,[], ...
%                 Meta_Data.PROCESS.nfft, ...
%                 Meta_Data.PROCESS.Fs_epsi,'psd');
% 
% % Create path, add shear Sv numbers, add CTD calibration,
% % save CTD calibration files in ctd and Shear calibration files in epsi
% Meta_Data=mod_define_meta_data_log(Meta_Data);
% 
% L1path = Meta_Data.L1path;
% clear Meta_Data
load(fullfile(path_mission,'L1','Meta_Data'))
