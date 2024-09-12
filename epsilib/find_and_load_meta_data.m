function [obj] = find_and_load_som_setup(obj,data_path)
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
dir_has_log = dir(fullfile(data_path,'Log*.csv'));
dir_has_config = dir(fullfile(data_path,'*config*'));

if ~isempty(dir_has_log) %if there is a log file...
    
    try
        obj.Meta_Data = create_metadata_from_deployment_log_v2(dir_has_log.name);
        obj.Meta_Data.AFE=obj.Meta_Data.epsi;
    catch err
        
        for j = 1:length(err.stack)
            disp([num2str(j) ' ' err.stack(j).name ' ' num2str(err.stack(j).line)]);
        end
        error('Failed to find config data (1)')
    end
    
elseif ~isempty(dir_has_config) %if there is a config file...
    
    try
        setup=mod_som_read_setup_from_config(dir_has_config.name);
    catch err
        for j = 1:length(err.stack)
            disp([num2str(j) ' ' err.stack(j).name ' ' num2str(err.stack(j).line)]);
        end
        error('Failed to find config data (2)')
        
    end
    % Fill Meta Data from setup data
    try
        obj.Meta_Data = epsiSetup_fill_meta_data(obj.Meta_Data,setup);
        
        fprintf('Meta_Data.paths.process_library is %s \n',obj.Meta_Data.paths.process_library);
        fprintf('Meta_Data.paths.data is %s \n',obj.Meta_Data.paths.data);
    catch err
        
        for j = 1:length(err.stack)
            disp([num2str(j) ' ' err.stack(j).name ' ' num2str(err.stack(j).line)]);
        end
        error('fill_meta_data failed (2)')
    end
    
else %if there is no log file or config file, look for config data inside the raw files
    
    % Try loading the most recent file in the raw directory
    try
        setupfile=dir(fullfile(obj.Meta_Data.paths.raw_data,...
            ['*' obj.Meta_Data.rawfileSuffix]));
        setup=mod_som_read_setup_from_raw(fullfile(setupfile(end).folder,setupfile(end).name));
    catch err
        for j = 1:length(err.stack)
            disp([num2str(j) ' ' err.stack(j).name ' ' num2str(err.stack(j).line)]);
        end
        warning(['Failed to read config data (3) - '...
            'this is often because mod_som_read_setup_from raw does not '...
            'have the correct offsets and lengths. When changes are made on '...
            'the hardware side, they have to be made here too.'])
        
        % Copy bench_config to data directory
        sourceFile = fullfile(obj.Meta_Data.paths.process_library,'config_files','bench_config');
        destinationDir = obj.Meta_Data.paths.data;
        copyfile(sourceFile, destinationDir);
        
        % Read setup data from bench_config
        dir_has_config = dir(fullfile(data_path,'*config*'));
        setup=mod_som_read_setup_from_config(dir_has_config.name);

        fprintf('Copied bench_config to data directory. \n')
    end

    % Fill Meta Data from setup data
    try
        obj.Meta_Data = epsiSetup_fill_meta_data(obj.Meta_Data,setup);
        
        fprintf('Meta_Data.paths.process_library is %s \n',obj.Meta_Data.paths.process_library);
        fprintf('Meta_Data.paths.data is %s \n',obj.Meta_Data.paths.data);
    catch err
        for j = 1:length(err.stack)
            disp([num2str(j) ' ' err.stack(j).name ' ' num2str(err.stack(j).line)]);
        end
        error('fill_meta_data failed (3)')
    end
    
end


