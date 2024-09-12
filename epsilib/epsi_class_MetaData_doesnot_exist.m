function [repeat,obj]=epsi_class_MetaData_doesnot_exist(data_path,obj)

%fprintf('Initializing epsi_class and creating new Meta_Data \n')
repeat = 0;

%%
% Define data path
obj.Meta_Data.paths.data     = data_path;
obj.Meta_Data.paths.raw_data = fullfile(data_path,'raw');
obj.Meta_Data.paths.mat_data = fullfile(data_path,'mat');
obj.Meta_Data.paths.profiles = fullfile(data_path,'profles');
obj.Meta_Data.paths.figures  = fullfile(data_path,'figs');
%%
% Find the epsi library and add it as process path
spltpath=strsplit(path,':');
epsilib_path=spltpath{~cellfun(@isempty, ...
    cellfun(@(x) ...
    strfind(x,'epsilib'),spltpath, ...
    'UniformOutput',false))};

obj.Meta_Data.paths.process_library=fileparts(epsilib_path);

addpath(genpath(obj.Meta_Data.paths.process_library));

rmpath(genpath(fullfile(obj.Meta_Data.paths.process_library,'archived_scripts')))

%%
% Add calibrations path
obj.Meta_Data.paths.calibration = fullfile(obj.Meta_Data.paths.process_library,'CALIBRATION','ELECTRONICS');
%%
% Read PROCESS Meta_Data from text file -
% if one is not specified, use the default
if isempty(obj.Meta_Data.paths.process_library)
    if isclassfield(obj.Meta_Data,'PROCESS')

        if ~isclassfield(obj.Meta_Data.PROCESS,'filename')
            Meta_Data_process_file = fullfile(obj.Meta_Data.paths.process_library,'Meta_Data_Process','Meta_Data_Process.txt');
        else
            % Use the .txt file save in Meta_Data, but find it
            % in the current user's path
            [~,fname,fsuffix] = fileparts(obj.Meta_Data.PROCESS.filename);
            Meta_Data_process_file = fullfile(obj.Meta_Data.paths.process_library,...
                'Meta_Data_Process',[fname,fsuffix]);
        end
    else
        Meta_Data_process_file = fullfile(obj.Meta_Data.paths.process_library,'Meta_Data_Process','Meta_Data_Process.txt');
    end

else
    % ALB Nicole has the following line I do not quite understand it. 
    % Meta_Data_process_file = Meta_Data_process_file;
    Meta_Data_process_file = fullfile(obj.Meta_Data.paths.process_library,'Meta_Data_Process','Meta_Data_Process.txt');
end


obj.f_read_MetaProcess(Meta_Data_process_file);

obj.Meta_Data = epsiSetup_set_epsi_paths(obj.Meta_Data);
obj.Meta_Data = epsiSetup_get_raw_suffix(obj.Meta_Data);


% NC TO DO:
% ADD A 4TH OPTION TO READ SETUP AND METADATA FROM YAML FILE
%
% There are three cases for getting configuration data
% and Meta_Data
%   1) Binary (?) config info is in the first raw file (done by
%   typing settings.stream during data aquistion,
%   header is $SOM3). Extra Meta_Data for processing is
%   in a Meta_Data_Process file.
%   2) Binary (?) config info is in a file called
%   *config*. Extra Meta_Data for processing is
%   in a Meta_Data_Process file.
%   3) All Meta_Data is in a csv file called Log_*.csv

% Is there a log csv file? Is there a config file? Or
% is config data inside the raw data files?
% ALB 2024/08/23 We are trying to get rid of log and bench_config
% ALB RIght I am commenting the following and grab everything from the
% .modraw everytime we read a new one.

% dir_has_log = dir(fullfile(data_path,'Log*.csv'));
% dir_has_config = dir(fullfile(data_path,'*config*'));

% if ~isempty(dir_has_log) %if there is a log file...
% 
%     try
%         obj.Meta_Data = create_metadata_from_deployment_log_v2(dir_has_log.name);
%         obj.Meta_Data.AFE=obj.Meta_Data.epsi;
%     catch err
% 
%         for j = 1:length(err.stack)
%             disp([num2str(j) ' ' err.stack(j).name ' ' num2str(err.stack(j).line)]);
%         end
%         warning('Failed to find config data (1)')
%         return
%     end
% 
% elseif ~isempty(dir_has_config) %if there is a config file...
%     try
%         % Sometimes there is a hidden file
%         % '._bench_config'. Get rid of that if it's there.
%         % You can just use dir_has_config(end).name and
%         % that will always take the last file in that
%         % directory that matches '*config*'
%         setup=mod_som_read_setup_from_config(dir_has_config(end).name);
%     catch err
%         for j = 1:length(err.stack)
%             disp([num2str(j) ' ' err.stack(j).name ' ' num2str(err.stack(j).line)]);
%         end
%         error('Failed to find config data (2)')
% 
%     end

    % % ALB this could be removed when reading Meta_Data
    % % from modraw (if SBE cal coef and epsi cal coef are printed there)
    % % NC added on NORSE 2023 - add serial numbers
    % % defined in Meta_Data_Process file
    % try
    %     setup.S49.sn = obj.Meta_Data.PROCESS.ctd_sn;
    %     setup.EFE.sensors{1}.sn = obj.Meta_Data.PROCESS.t1_sn;
    %     setup.EFE.sensors{2}.sn = obj.Meta_Data.PROCESS.t2_sn;
    %     setup.EFE.sensors{3}.sn = obj.Meta_Data.PROCESS.s1_sn;
    %     setup.EFE.sensors{4}.sn = obj.Meta_Data.PROCESS.s2_sn;
    % catch
    % end

    % % Fill Meta Data from setup data
    % try
    %     obj.Meta_Data = epsiSetup_fill_meta_data(obj.Meta_Data,setup);
    % 
    %     %fprintf('Meta_Data.paths.process_library is %s \n',obj.Meta_Data.paths.process_library);
    %     %fprintf('Meta_Data.paths.data is %s \n',obj.Meta_Data.paths.data);
    % catch err
    % 
    %     for j = 1:length(err.stack)
    %         disp([num2str(j) ' ' err.stack(j).name ' ' num2str(err.stack(j).line)]);
    %     end
    %     warning('fill_meta_data failed (2)')
    %     return
    % end

% else %if there is no log file or config file, look for config data inside the raw files
    % TODO 10/7/21 - Loop through more that just the
    % first file to look for $SOM3

    % ALB also getting of this as we should read do this every time we open
    % a modraw
    % try
    %     setupfile=dir(fullfile(obj.Meta_Data.paths.raw_data,...
    %         ['*' obj.Meta_Data.rawfileSuffix]));
    %     setup=mod_som_read_setup_from_raw(fullfile(setupfile(1).folder,setupfile(1).name));
    % catch err
    %     for j = 1:length(err.stack)
    %         disp([num2str(j) ' ' err.stack(j).name ' ' num2str(err.stack(j).line)]);
    %     end
    %     warning(['Failed to read config data (3) - '...
    %         'this is often because mod_som_read_setup_from raw does not '...
    %         'have the correct offsets and lengths. When changes are made on '...
    %         'the hardware side, they have to be made here too.'])
    % 
    %     % Copy bench_config to data directory
    %     sourceFile = fullfile(obj.Meta_Data.paths.process_library,'config_files','bench_config');
    %     destinationDir = obj.Meta_Data.paths.data;
    %     copyfile(sourceFile, destinationDir);
    % 
    %     % Read setup data from bench_config
    %     dir_has_config = dir(fullfile(data_path,'*config*'));
    %     setup=mod_som_read_setup_from_config(dir_has_config.name);
    % 
    %     fprintf('Copied bench_config to data directory. \n')
    % end

    % ALB this could be removed when reading Meta_Data
    % from modraw (if SBE cal coef and epsi cal coef are printed there)
    % NC added on NORSE 2023 - add serial numbers
    % defined in Meta_Data_Process file
    % try
    %     setup.S49.sn = obj.Meta_Data.PROCESS.ctd_sn;
    %     setup.EFE.sensors{1}.sn = obj.Meta_Data.PROCESS.t1_sn;
    %     setup.EFE.sensors{2}.sn = obj.Meta_Data.PROCESS.t2_sn;
    %     setup.EFE.sensors{3}.sn = obj.Meta_Data.PROCESS.s1_sn;
    %     setup.EFE.sensors{4}.sn = obj.Meta_Data.PROCESS.s2_sn;
    % catch
    % end

    % % Fill Meta Data from setup data
    % try
    %     obj.Meta_Data = epsiSetup_fill_meta_data(obj.Meta_Data,setup);
    % 
    %     %fprintf('Meta_Data.paths.process_library is %s \n',obj.Meta_Data.paths.process_library);
    %     %fprintf('Meta_Data.paths.data is %s \n',obj.Meta_Data.paths.data);
    % catch err
    %     for j = 1:length(err.stack)
    %         disp([num2str(j) ' ' err.stack(j).name ' ' num2str(err.stack(j).line)]);
    %     end
    %     warning('fill_meta_data failed (3)')
    %     return
    % end

% end

% % Set epsi paths and define suffix for raw files
% % Read PROCESS Meta_Data from default text file
% obj.f_read_MetaProcess(Meta_Data_process_file);

% % ALB this could be removed when reading Meta_Data
% % from modraw (if SBE cal coef and epsi cal coef are printed there)
% % NC added on NORSE 2023 - add serial numbers
% % defined in Meta_Data_Process file
% try
%     setup.S49.sn = obj.Meta_Data.PROCESS.ctd_sn;
%     setup.EFE.sensors{1}.sn = obj.Meta_Data.PROCESS.t1_sn;
%     setup.EFE.sensors{2}.sn = obj.Meta_Data.PROCESS.t2_sn;
%     setup.EFE.sensors{3}.sn = obj.Meta_Data.PROCESS.s1_sn;
%     setup.EFE.sensors{4}.sn = obj.Meta_Data.PROCESS.s2_sn;
% catch
% end

% obj.Meta_Data = epsiSetup_set_epsi_paths(obj.Meta_Data);
% obj.Meta_Data = epsiSetup_get_raw_suffix(obj.Meta_Data);
Meta_Data = obj.Meta_Data;
save(fullfile(obj.Meta_Data.paths.data,'Meta_Data'),'Meta_Data');

