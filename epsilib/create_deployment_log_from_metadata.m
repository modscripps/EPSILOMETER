% In order to make epsi_class backwards compatible with older datasets that
% don't use a config file, I want to make sure each of these old
% deployments has a csv log file. When epsi_class initializes, it will look
% for either a csv file, or a config file, or the config data inside one of
% the raw files.
% - Use this function to create log files for old deployments that already
% have Meta_Data.mat.
% - Other Meta_Data.PROCESS information will go into a
% Meta_Data_Process.txt file

function [] = create_deployment_log_from_metadata(input_file_path,output_file_name)

% Load Meta_Data
load(input_file_path);

% Define output file location and name
[dir1,~]=fileparts(input_file_path);
[dir2,~]=fileparts(dir1);
output_file = fullfile(dir2,...
    output_file_name);

parameters{1} = 'Mission name';
try
    values{1} = Meta_Data.mission;
catch
    values{1} = ' ';
end

parameters{2} = 'Path mission';
try
    values{2} = Meta_Data.path_mission;
catch
    values{2} = ' ';
end

parameters{3} = 'Vehicle';
try
    values{3} = Meta_Data.vehicle;
catch
    values{3} = ' ';
end

parameters{4} = 'Vehicle name';
try
    values{4} = Meta_Data.vehicle_name;
catch
    values{4} = ' ';
end

parameters{5} = 'Deployment name';
try
    values{5} = Meta_Data.deployment;
catch
    values{5} = ' ';
end

parameters{6} = 'NFFT';
try
    values{6} = num2str(Meta_Data.PROCESS.nfft);
catch
    values{6} = ' ';
end

parameters{7} = 'Fs epsi';
try
    values{7} = num2str(Meta_Data.PROCESS.Fs_epsi);
catch
    values{7} = ' ';
end

parameters{8} = 'Shear probe 1';
try
    values{8} = Meta_Data.epsi.s1.SN;
catch
    values{8} = ' ';
end

parameters{9} = 'Shear probe 2';
try
    values{9} = Meta_Data.epsi.s2.SN;
catch
    values{9} = ' ';
end
parameters{10} = 'FPO7 1';
try
    values{10} = Meta_Data.epsi.t1.SN;
catch
    values{10} = ' ';
end

parameters{11} = 'FPO7 2';
try
    values{11} = Meta_Data.epsi.t2.SN;
catch
    values{11} = ' ';
end

parameters{12} = 'CTD name';
try
    values{12} = Meta_Data.CTD.name;
catch
    values{12} = ' ';
end

parameters{13} = 'CTD SN';
try
    values{13} = Meta_Data.CTD.SN;
catch
    values{13} = ' ';
end

parameters{14} = 'Channels';
try
    values{14} = string(Meta_Data.PROCESS.channels);
catch
    values{14} = ' ';
end

parameters{15} = 'Recording mode';
try
    values{15} = Meta_Data.PROCESS.recording_mode;
catch
    values{15} = ' ';
end

parameters{16} = 'Starttime';
try
    values{16} = datestr(Meta_Data.starttime,'dd-mmm-yyyy HH:MM:SS');
catch
    values{16} = '00-Jan-0000';
end

parameters{17} = 'MADRE rev';
try
    values{17} = Meta_Data.MADRE.rev;
catch
    values{17} = ' ';
end

parameters{18} = 'MADRE SN';
try
    values{18} = Meta_Data.MADRE.SN;
catch
    values{18} = ' ';
end
parameters{19} = 'MAP rev';
try
    values{19} = Meta_Data.MAP.rev;
catch
    values{19} = ' ';
end

parameters{20} = 'MAP SN';
try
    values{20} = Meta_Data.MAP.SN;
catch
    values{20} = ' ';
end

parameters{21} = 'MAP temperature';
try
    values{21} = Meta_Data.MAP.temperature;
catch
    values{21} = ' ';
end

parameters{22} = 'MAP shear';
try
    values{22} = Meta_Data.MAP.shear;
catch
    values{22} = ' ';
end

parameters{23} = 'Firmware version';
try
    values{23} = Meta_Data.Firmware.version;
catch
    values{23} = ' ';
end

parameters{24} = 'ADC shear';
try
    values{24} = Meta_Data.Firmware.ADCshear;
catch
    values{24} = ' ';
end

parameters{25} = 'ADC FPO7';
try
    values{25} = Meta_Data.Firmware.ADC_FPO7;
catch
    values{25} = ' ';
end

parameters{26} = 'ADC accelerometer';
try
    values{26} = Meta_Data.Firmware.ADC_accellerometer;
catch
    values{26} = ' ';
end

parameters{27} = 'ADC conductivity';
try
    values{27} = Meta_Data.Firmware.ADC_cond;
catch
    values{27} = ' ';
end

parameters{28} = 'ADC sampling filter';
try
    values{28} = Meta_Data.Firmware.ADCfilter;
catch
    values{28} = ' ';
end

% Save Meta_Data in a table, MD
MD = table(parameters(:),values(:));

% Write data to output_file
writetable(MD,output_file,'delimiter',',')

fid = fopen(output_file,'w');
fprintf(fid,'%s,%s,\n','parameter','value')
for ii=1:length(parameters)
    if ii~=14
        fprintf(fid,'%s,%s,\n',parameters{ii},values{ii});
    elseif ii==14
        formatString = repmat('%s,',1,length(values{ii}));
        fprintf(fid,['%s,' formatString '\n'],parameters{ii},values{ii});
    end
end
fclose(fid);

