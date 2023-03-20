classdef epsi_class < handle
    %
    % epsi_class is a class of epsi data, processing functions, and
    % plotting functions
    %
    % Usage:
    % To create an epsi_class structure, you must first be within a
    % directory containing 1) an epsi config file, and 2) a directory
    % called 'raw' containing raw files which all have the same file
    % suffix (ex. .ascii, .raw, _raw). Make sure the suffix is one of the
    % options in epsiSetup_get_raw_suffix.m
    %   >> obj = epsi_class(lastOrAll)
    %
    % OPTIONAL INPUTS
    %   data_path = path to data directory (the parent directory that
    %               contains raw, mat, etc)
    %   Meta_Data_Process_file = path to .txt file with Meta_Data.PROCESS
    %                            values

    properties
        Meta_Data %read from config file
        filedir
        plot_properties
        epsi
        ctd
        alt
        vnav
        gps
    end
    methods
        function obj=epsi_class(varargin)
            if nargin<1
                data_path = pwd;
                Meta_Data_process_file = []; %Default Meta_Data process file not specified
            elseif nargin==1 %varargin is either data directory or Meta_Data_Process_file
                if exist(varargin{1},'dir') %if varargin is a directory
                    data_path = varargin{1};
                    Meta_Data_process_file = [];
                elseif exist(varargin{1},'file') %if it's a file
                    Meta_Data_process_file = varargin{1};
                    data_path = pwd;
                end

            elseif nargin==2 %varagin is data directory and Meta_Data_Process file
                if exist(varargin{1},'dir') %if the first input is a directory
                    data_path = varargin{1};
                    Meta_Data_process_file = varargin{2};
                elseif exist(varargin{1},'file') %if the first input is a file
                    Meta_Data_process_file = varargin{1};
                    data_path = varargin{2};
                end
            end

            %             % removed archive+path
            %             spltpath=strsplit(path,':');
            %             try
            %                 archived_path=spltpath{~cellfun(@isempty, ...
            %                     cellfun(@(x) ...
            %                     strfind(x,'archived_scripts'),spltpath, ...
            %                     'UniformOutput',false))};
            %                 if ~isempty(archived_path)
            %                     rmpath(genpath(archived_path));
            %                 end
            %             catch
            %             end

            % Check to see if Meta_Data is already defined
            checkMD = dir(fullfile(data_path,'Meta_Data.mat'));

            repeat = 1; %Initialize repeat flag to use if Meta_Data path names were not made on this machine
            while repeat==1

                if ~isempty(checkMD) %Meta_Data already exists

                    fprintf('Initializing epsi_class with previously created Meta_Data \n')
                    load(fullfile(data_path,'Meta_Data'))
                    obj.Meta_Data = Meta_Data;

                    % Meta_Data includes the path to the epsi processing
                    % library, but if the data was processed on another
                    % machine, you won't have access to it.
                    % Find the epsi library on your machine and add it as process path
                    spltpath=strsplit(path,':');
                    epsilib_path=spltpath{~cellfun(@isempty, ...
                        cellfun(@(x) ...
                        strfind(x,'epsilib'),spltpath, ...
                        'UniformOutput',false))};
                    obj.Meta_Data.paths.process_library=fileparts(epsilib_path);
                    obj.Meta_Data.paths.calibration = fullfile(obj.Meta_Data.paths.process_library,'CALIBRATION','ELECTRONICS');

                    % Always redefine the data path as the current
                    % directory or the directory you input
                    obj.Meta_Data.paths.data=data_path;
                    obj.Meta_Data.paths.raw_data = fullfile(data_path,'raw');
                    obj.Meta_Data.paths.mat_data = fullfile(data_path,'mat');
                    obj.Meta_Data.paths.profiles = fullfile(data_path,'profles');
                    obj.Meta_Data.paths.figures = fullfile(data_path,'figs');

                    obj.Meta_Data = epsiSetup_get_raw_suffix(obj.Meta_Data);

                    % Check that Meta_Data has everything you need. If it's
                    % missing something, set repeat=0 and checkMD=[]
                    if  isdir(obj.Meta_Data.paths.process_library) && ...
                            isdir(obj.Meta_Data.paths.data) && ...
                            isdir(obj.Meta_Data.paths.calibration) && ...
                            isclassfield(obj.Meta_Data.paths,'raw_data') && ...
                            isclassfield(obj.Meta_Data.PROCESS,'nb_channels') && ...
                            isclassfield(obj.Meta_Data.PROCESS,'nfft')
                        Meta_Data = obj.Meta_Data;
%                         save(fullfile(obj.Meta_Data.paths.data,'Meta_Data'),'Meta_Data');

                        repeat = 0; %Stop repeating. Keep this Meta_Data.
                    else
                        checkMD = []; %Some parts are missing. Repeat to make new Meta_Data.
                    end

                elseif isempty(checkMD) %Meta_Data, epsi, and ctd .mat files do not already exist

                    fprintf('Initializing epsi_class and creating new Meta_Data \n')
                    repeat = 0;

                    % Define data path
                    obj.Meta_Data.paths.data     = data_path;
                    obj.Meta_Data.paths.raw_data = fullfile(data_path,'raw');
                    obj.Meta_Data.paths.mat_data = fullfile(data_path,'mat');
                    obj.Meta_Data.paths.profiles = fullfile(data_path,'profles');
                    obj.Meta_Data.paths.figures  = fullfile(data_path,'figs');

                    % Find the epsi library and add it as process path
                    spltpath=strsplit(path,':');
                    epsilib_path=spltpath{~cellfun(@isempty, ...
                        cellfun(@(x) ...
                        strfind(x,'epsilib'),spltpath, ...
                        'UniformOutput',false))};
                    obj.Meta_Data.paths.process_library=fileparts(epsilib_path);
                    addpath(genpath(obj.Meta_Data.paths.process_library));
                    rmpath(genpath(fullfile(obj.Meta_Data.paths.process_library,'archived_scripts')))

                    % Add calibrations path
                    obj.Meta_Data.paths.calibration = fullfile(obj.Meta_Data.paths.process_library,'CALIBRATION','ELECTRONICS');

                    % Read PROCESS Meta_Data from text file -
                    % if one is not specified, use the default
                    if isempty(Meta_Data_process_file)
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
                        Meta_Data_process_file = Meta_Data_process_file;
                    end


                    obj.f_read_MetaProcess(Meta_Data_process_file);

                    obj.Meta_Data = epsiSetup_set_epsi_paths(obj.Meta_Data);
                    obj.Meta_Data = epsiSetup_get_raw_suffix(obj.Meta_Data);

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
                        % TODO 10/7/21 - Loop through more that just the
                        % first file to look for $SOM3

                        try
                            setupfile=dir(fullfile(obj.Meta_Data.paths.raw_data,...
                                ['*' obj.Meta_Data.rawfileSuffix]));
                            setup=mod_som_read_setup_from_raw(fullfile(setupfile(1).folder,setupfile(1).name));
                        catch err
                            for j = 1:length(err.stack)
                                disp([num2str(j) ' ' err.stack(j).name ' ' num2str(err.stack(j).line)]);
                            end
                            error(['Failed to read config data (3) - '...
                                'this is often because mod_som_read_setup_from raw does not '...
                                'have the correct offsets and lengths. When changes are made on '...
                                'the hardware side, they have to be made here too.'])
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

                    % Set epsi paths and define suffix for raw files
                    % Read PROCESS Meta_Data from default text file
                    obj.f_read_MetaProcess(Meta_Data_process_file);
                    obj.Meta_Data = epsiSetup_set_epsi_paths(obj.Meta_Data);
                    obj.Meta_Data = epsiSetup_get_raw_suffix(obj.Meta_Data);
                    Meta_Data = obj.Meta_Data;
                    save(fullfile(obj.Meta_Data.paths.data,'Meta_Data'),'Meta_Data');

                end
            end

            % Always reset epsi paths to current directory. Otherwise you
            % could end up reprocessing things in an old directory you
            % tried to move away from (personal experiece)
            obj.Meta_Data = epsiSetup_set_epsi_paths(obj.Meta_Data);

            % Also, print all the directories just to be sure
            disp('... These are your directories:')
            fprintf('process path: %s \n',obj.Meta_Data.paths.process_library);
            fprintf('data path: %s \n',obj.Meta_Data.paths.data);
            fprintf('calibration path: %s \n',obj.Meta_Data.paths.calibration);
            fprintf('raw data path: %s \n',obj.Meta_Data.paths.raw_data);
            fprintf('mat data path: %s \n',obj.Meta_Data.paths.mat_data);
            fprintf('profiles path: %s \n',obj.Meta_Data.paths.profiles);

            % NC - check for s1 and s2 cal values. For newer deployments that have Meta_Data.AFE structure, if they're
            % 0, manually input all probe numbers.
            % NC 10/7/21 - Check for 'AFE' or 'epsi' strucutre in Meta_Data. Add
            % calibratation to the appropriate structure.
            if isclassfield(obj.Meta_Data,'AFE') && ~isclassfield(obj.Meta_Data,'epsi')
                field_name = 'AFE';
            elseif isclassfield(obj.Meta_Data,'epsi') && ~isclassfield(obj.Meta_Data,'AFE')
                field_name = 'epsi';
            elseif isclassfield(obj.Meta_Data,'epsi') && isclassfield(obj.Meta_Data,'AFE')
                field_name = 'epsi';
            else
                field_name = [];
            end

            if ~isempty(field_name)
                obj.Meta_Data = obj.f_getSNshear;
                obj.Meta_Data = obj.f_getSNtemp;
            end

            % Read PROCESS Meta_Data from default text file -
            % if one is not specified, use the default
            if isempty(Meta_Data_process_file)  && ~isclassfield(Meta_Data.PROCESS,'filename')
                Meta_Data_process_file = fullfile(obj.Meta_Data.paths.process_library,'Meta_Data_Process','Meta_Data_Process.txt');
            elseif isclassfield(Meta_Data.PROCESS,'filename')
                Meta_Data_process_file = obj.Meta_Data.PROCESS.filename;
            end
            obj.f_read_MetaProcess(Meta_Data_process_file);

            % Define filedir as path to raw data
            obj.filedir=obj.Meta_Data.paths.data;

            % Define plot properties
            obj.f_getPlotProperties;
        end
        function obj=f_read_MetaProcess(obj,filename)
            if nargin==1
                filename=fullfile(obj.Meta_Data.paths.process_library,'Meta_Data_Process',...
                    'Meta_Data_Process.txt');
            end
            try
                Meta_Data = epsiSetup_read_MetaProcess(obj.Meta_Data,filename);
            catch
                Meta_Data = load(filename);
                Meta_Data=Meta_Data.Meta_Data;
            end
            obj.Meta_Data = Meta_Data;
        end
        function obj=f_readData(obj,varargin)
            % optional arguments:
            %   'version'
            %       ex) MISOBOB - f_readData('version',0)
            %       ex) BLT2021 - f_readData('version',3)
            %       ex) BLT2022 - f_readData('version',4)
            %   'calc_micro'
            %       ex) f_readData('calc_micro') will calculate microstructure
            %       ex) f_readData will not calculate microstructrure
            %   'make_FCTD'
            %       ex) f_readData('make_FCTD','path/to/fctd_mat_files')
            %                       will make FCTD-style .mat files in the
            %                       path specified in the second argument
            %       ex) f_readData('make_FCTD') will make FCTD-style .mat
            %                       files in a directory called FCTDmat at
            %                       the same level as .mat and .raw directories
            %       ex) f_readData will not not make FCTD-style .mat files


            % Set defaults
            version_number = 4;
            calc_micro = 0;
            make_FCTD = 0;
            fctd_mat_dir = '';

            argsNameToCheck = {'calc_micro',...     %1
                'version',...        %2
                'make_FCTD'};        %3

            index = 1; %Initialize index of argsNameToCheck
            % Number of items remaining (this is the number of argsNameToCheck minus
            % the number of extra parameters that go with the arguments. For example,
            % 'version' expects another parameter that follows it: 'version', 3.
            % Similarly, 'fileStr' expects a string after it: 'fileStr',
            % 'EPSI_22_04_12*'
            n_items = nargin-2;

            while (n_items > 0)
                argsMatch = strcmpi(varargin{index},argsNameToCheck);
                i = find(argsMatch,1);
                if isempty(i)
                    error('MATLAB:epsiProcess_convert_new_raw_to_mat:wrongOption','Incorrect option specified: %s', varargin{index});
                end

                switch i
                    case 1 %calc_micro
                        % Find the index of varargin that = 'calc_micro'
                        % and set calc_micro to 1
                        idxFlag = find(cell2mat(cellfun(@(C) ~isempty(strfind(C,'calc_micro')),varargin,'uniformoutput',0)));
                        calc_micro = 1;
                        index = index+1;
                        n_items = n_items-1;
                    case 2 %version
                        % Find the index of varargin that = 'version'. The following
                        % index contains the version number
                        idxFlag = find(cell2mat(cellfun(@(C) ~isempty(strfind(C,'version')),varargin,'uniformoutput',0)));
                        version_number = varargin{idxFlag+1};
                        index = index+2; %+2 because the following varargin will be the version number
                        n_items = n_items-2;
                    case 3 %make_FCTD
                        % Find the index of varargin that = 'make_FCTD' and
                        % set make_FCTD to 1
                        idxFlag = find(cell2mat(cellfun(@(C) ~isempty(strfind(C,'make_FCTD')),varargin,'uniformoutput',0)));
                        make_FCTD = 1;
                        index = index+1;
                        n_items = n_items-1;
                end
            end

            % Copy raw files from datapath to RAWpath
            %list_rawfile = dir(fullfile(obj.Meta_Data.paths.data,['*' obj.Meta_Data.rawfileSuffix]));
            list_rawfile = dir(fullfile(obj.Meta_Data.paths.raw_data,['*' obj.Meta_Data.rawfileSuffix]));

            % Some files are called '*data*' and if you look for them, you
            % might  also grab Meta_Data. Get rid of it.
            matchCell = strfind({list_rawfile.name}, 'Meta_Data');
            matchArray = cellfun(@(C) ~isempty(C),matchCell);
            list_rawfile = list_rawfile(~matchArray);

            if isempty(list_rawfile)
                warning(['There are no *' obj.Meta_Data.rawfileSuffix ' * raw files in ' obj.Meta_Data.paths.raw_data])
            else
                %                 for f=1:length(list_rawfile)
                %                     copyfile(fullfile(list_rawfile(f).folder, ...
                %                         list_rawfile(f).name),  ...
                %                         obj.Meta_Data.paths.raw_data,'f'); %NC set mode to 'f' to copy file
                %                 end
                % Convert raw to mat
                dirs.raw_incoming = obj.Meta_Data.paths.raw_data;
                dirs.mat = obj.Meta_Data.paths.mat_data;
                if make_FCTD
                    epsiProcess_convert_new_raw_to_mat(dirs,obj.Meta_Data,...
                        'noSync',...
                        'version',version_number,...
                        'calc_micro',calc_micro,...
                        'make_FCTD',fctd_mat_dir);
                elseif ~make_FCTD
                    epsiProcess_convert_new_raw_to_mat(dirs,obj.Meta_Data,...
                        'noSync',...
                        'version',version_number,...
                        'calc_micro',calc_micro);
                end
            end

        end
        function obj=f_getLastData(obj)
            obj.f_readData;
            obj.epsi = obj.f_getLastEpsi();
            obj.ctd = obj.f_getLastCtd();
            obj.alt = obj.f_getLastAlt();
        end
        function obj=f_getLastEpsi(obj)
            load(fullfile(obj.Meta_Data.paths.mat_data,'TimeIndex'));
            [~,idxLast] = max(TimeIndex.timeEnd);
            data = load(fullfile(obj.Meta_Data.paths.mat_data,TimeIndex.filenames{idxLast}),'epsi');

            if isstruct(data.epsi)
                obj=data.epsi;
            else
                obj=[];
            end
        end
        function obj=f_getLastCtd(obj)
            load(fullfile(obj.Meta_Data.paths.mat_data,'TimeIndex'));
            [~,idxLast] = max(TimeIndex.timeEnd);
            data = load(fullfile(obj.Meta_Data.paths.mat_data,[TimeIndex.filenames{idxLast} '.mat']),'ctd');

            if isstruct(data.ctd)
                obj=data.ctd;
            else
                obj=[];
            end
        end
        function obj=f_getLastAlt(obj)
            load(fullfile(obj.Meta_Data.paths.mat_data,'TimeIndex'));
            [~,idxLast] = max(TimeIndex.timeEnd);
            data = load(fullfile(obj.Meta_Data.paths.mat_data,TimeIndex.filenames{idxLast}),'alt');

            if isstruct(data.alt)
                obj=data.alt;
            else
                obj=[];
            end
        end


        function obj=f_getFileData(obj,fileNumOrName)
            % fileNumOrName is an array of file indices of the list of raw
            % files Or the name of the file
            obj.f_readData;

            load(fullfile(obj.Meta_Data.paths.mat_data,'TimeIndex'));

            if isnumeric(fileNumOrName)
                if length(fileNumOrName)==1 %If you chose one file, as an index number
                    idx = fileNumOrName;
                    fileList{1} = fullfile(obj.Meta_Data.paths.mat_data,TimeIndex.filenames{idx});
                else %If you chose more than one file as a list of indices
                    for iF=1:length(fileNumOrName)
                        fileList{iF} = fullfile(obj.Meta_Data.paths.mat_data,TimeIndex.filenames{iF});
                    end
                end
            elseif ischar(fileNumOrName) %If you chose one file as a character string
                try
                    fileList{1} = fullfile(obj.Meta_Data.paths.mat_data,fileNumOrName);
                    load(fileList{1})
                catch
                    fileList{1} = fileNumOrName;
                    load(fileList{1})
                end
            end

            data = epsiSetup_make_empty_structure();
            for iF=1:length(fileList)-1

                fileName = fileList{iF+1};

                data2.epsi = obj.f_getFileEpsi(fileName);
                data2.ctd = obj.f_getFileCtd(fileName);
                data2.alt = obj.f_getFileAlt(fileName);

                data.epsi = epsiProcess_merge_mat_files(data.epsi,data2.epsi);
                data.ctd  = epsiProcess_merge_mat_files(data.ctd,data2.ctd);
                data.alt  = epsiProcess_merge_mat_files(data.alt,data2.alt);
            end

            obj.epsi = data.epsi;
            obj.ctd = data.ctd;
            obj.alt = data.alt;
        end
        function obj=f_getFileEpsi(obj,fileName)
            data = load(fileName,'epsi');
            if isstruct(data.epsi)
                obj=data.epsi;
            else
                obj=[];
            end
        end
        function obj=f_getFileCtd(obj,fileName)
            data = load(fileName,'ctd');
            if isstruct(data.ctd)
                obj=data.ctd;
            else
                obj=[];
            end
        end
        function obj=f_getFileAlt(obj,fileName)
            data = load(fileName,'alt');
            if isstruct(data.alt)
                obj=data.alt;
            else
                obj=[];
            end
        end
        function var_timeseries = f_get_var_timeseries(obj,var_name)
            var_timeseries = epsiProcess_get_var_timeseries(obj,var_name);
        end
        function var_timeseries = f_plot_var_timeseries(obj,var_name);
            var_timeseries = epsiPlot_var_timeseries(obj,var_name);
        end
        function obj=f_getPlotProperties(obj)
            % Set default plot properties (fonts, colors, sizes)
            obj.plot_properties = epsiSetup_set_plot_properties;
        end
        function obj=f_getSNshear(obj)
            % Input shear probe serial numbers into Meta_Data and get
            % latest calibration values (Sv)
            %
            % USAGE
            %   obj.Meta_Data = f_getSNshear(obj)

            %             Meta_Data  = obj.Meta_Data;
            obj = epsiSetup_get_SN_shear(obj.Meta_Data);
        end
        function obj=f_getSNtemp(obj)
            % Input temperature probe serial numbers into Meta_Data and get
            % latest calibration values (Sv)
            %
            % USAGE
            %   obj.Meta_Data = f_getSNtemp(obj)

            %             Meta_Data  = obj.Meta_Data;
            obj = epsiSetup_get_SN_temp(obj.Meta_Data);
        end
        function f_plotCtd(obj)
            figure
            ax(1)=subplot(311);
            ax(2)=subplot(312);
            ax(3)=subplot(313);

            plot(ax(1),obj.ctd.time_s,obj.ctd.P,'k','LineWidth',obj.plot_properties.LineWidth)
            plot(ax(2),obj.ctd.time_s,obj.ctd.T,'r','LineWidth',obj.plot_properties.LineWidth)
            plot(ax(3),obj.ctd.time_s,obj.ctd.S,'b','LineWidth',obj.plot_properties.LineWidth)
            ax(1).XTickLabel='';
            ax(2).XTickLabel='';
            title(ax(1),sprintf('%s-%s-%s',strrep(obj.Meta_Data.mission,'_','\_'),...
                strrep(obj.Meta_Data.vehicle_name,'_','\_'),...
                strrep(obj.Meta_Data.deployment,'_','\_')));
            ylabel(ax(1),'P [dbar]','Fontsize',obj.plot_properties.FontSize,'FontName',obj.plot_properties.FontName)
            ylabel(ax(2),'T [˚C]','Fontsize',obj.plot_properties.FontSize,'FontName',obj.plot_properties.FontName)
            ylabel(ax(3),'S [psu]','Fontsize',obj.plot_properties.FontSize,'FontName',obj.plot_properties.FontName)
            ax(3).XLabel.String = 'seconds';
            set(ax(1),'Fontsize',obj.plot_properties.FontSize,'FontName',obj.plot_properties.FontName,'YDir','reverse')
            set(ax(2),'Fontsize',obj.plot_properties.FontSize,'FontName',obj.plot_properties.FontName)
            set(ax(3),'Fontsize',obj.plot_properties.FontSize,'FontName',obj.plot_properties.FontName)


        end
        function f_plot_epsi_detrended(obj)
            plot_epsi_detrended(obj);
        end
        function f_plotEpsi(obj)
            plot_epsi(obj);
        end
        function f_plotAlti(obj)
            figure
            if isclassfield(obj,'ctd') %NC equivalent of is field for class structure
                ax(1) = axes('Position',[0.14    0.15    0.7    0.7]);
                plot(obj.alt.alttime,obj.alt.dst,'.r');
                ax(1).YColor = 'r';
                ax(1).YLabel.String = 'altimeter';
                datetick(ax(1),'x','mm/dd HH:MM','keeplimits')
                ax(1).YAxisLocation = 'right';

                ax(2) = axes('Position',ax(1).Position);
                plot(obj.ctd.time_s,obj.ctd.P,'Color',obj.plot_properties.Colors.P);
                ax(2).YColor = obj.plot_properties.Colors.P;
                ax(2).Color = 'none';
                ax(2).YLabel.String = 'ctd pressure';
                ax(2).YDir = 'reverse';
                ax(2).YLim(1) = -2;
                ax(2).XTickLabel = '';

                lp = linkprop([ax(:)],'xlim');
            else
                plot(obj.alt.alttime,obj.alt.dst,'.r');
            end
        end
        function f_plotPressureTimeseries(obj)
            fileName = fullfile(obj.Meta_Data.paths.mat_data,'PressureTimeseries.mat');
            load(fileName);
            figure

            plot(PressureTimeseries.dnum,PressureTimeseries.P,'k');
            hold on
            s = plot(PressureTimeseries.dnum(PressureTimeseries.startprof),PressureTimeseries.P(PressureTimeseries.startprof),'og');
            e = plot(PressureTimeseries.dnum(PressureTimeseries.endprof),PressureTimeseries.P(PressureTimeseries.endprof),'sr');
            set(gca,'ydir','reverse')
            legend([s,e],{'start prof','end prof'})
            datetick(gca,'x')
            set(gca,'xticklabelrotation',45)
            title([datestr(nanmin(PressureTimeseries.dnum)) ' - ' ...
                datestr(nanmax(PressureTimeseries.dnum))]);
            grid on
        end
        function f_plotFallSpeed(obj)
            figure
            plot(obj.ctd.time_s,obj.ctd.dPdt,'.')
            hold on
            plot(obj.ctd.time_s,movmean(obj.ctd.dPdt,20),'.')
            xlabel('time_s (seconds)')
            ylabel('dPdt')
        end
        function f_plotEpsiCtd(obj,saveFig)
            if nargin<2
                saveFig=0;
            end
            epsiPlot_epsi_ctd_alt_timeseries(obj,saveFig)
        end


        function  [P11,f,noise,ax]=f_plot_spectra_at_tMid(obj,tmid,tscan,nSec,makeFig,saveFig,replaceData,ax)
            % Plots 30-sec timeseries from all channels and spectra from
            % user-defined tscan
            %
            % USAGE
            %   [P11,f] = obj.f_calibrateEpsi(tmid,tscan,makeFig,saveFig)
            %   Equivalent to:
            %   [P11,f] = mod_som_calibrate_epsi_tMid(obj,tmid,tscan,makeFig);
            %
            % INPUTS
            %   tmid = midpoint of scan (seconds)
            %   tscan = length of scan (seconds)
            %   nSec      - length of timeseries to plot
            %   makeFig = (OPTIONAL, flag to plot figure [0/1], default=1)
            %   saveFig = (OPTIONAL, flag to save figure [0/1], default=1)
            %
            % OUTPUTS (OPTIONAL)
            %   P11 = structure of frequency spectra for each channel
            %   f = frequency array
            %   noise = structure of shear and fpo7 noise coefficients
            if nargin<7
                ax = [];
                if nargin<6
                    replaceData=0;
                    if nargin<5
                        makeFig = 1;
                        saveFig = 0;
                    end
                    if nargin==5
                        saveFig = 1;
                    end
                end
            end

            [P11,f,noise,ax] = epsiPlot_spectra_at_tMid(obj,tmid,tscan,nSec,makeFig,saveFig,replaceData,ax);
        end


        function  [P11,f,noise]=f_calibrateEpsi(obj,tmid,tscan,makeFig,saveFig,replaceData,ax)
            % Plots 30-sec timeseries from all channels and spectra from
            % user-defined tscan
            %
            % USAGE
            %   [P11,f] = obj.f_calibrateEpsi(tmid,tscan,makeFig,saveFig)
            %   Equivalent to:
            %   [P11,f] = mod_som_calibrate_epsi_tMid(obj,tmid,tscan,makeFig);
            %
            % INPUTS
            %   tmid = midpoint of scan (seconds),
            %   default:  tmid= time_s(end) - 10.
            %   tscan = length of scan (seconds)
            %   makeFig = (OPTIONAL, flag to plot figure [0/1], default=1)
            %   saveFig = (OPTIONAL, flag to save figure [0/1], default=1)
            %
            % OUTPUTS (OPTIONAL)
            %   P11 = structure of frequency spectra for each channel
            %   f = frequency array
            %   noise = structure of shear and fpo7 noise coefficients
            if nargin<6
                ax = [];
                if nargin<5
                    replaceData=0;
                    if nargin<4
                        makeFig = 1;
                        saveFig = 0;
                    end
                    if nargin==4
                        saveFig = 1;
                    end
                end
            end
            [P11,f,noise,ax] = mod_som_calibrate_epsi_tMid(obj,tmid,tscan,makeFig,saveFig,replaceData,ax);
        end
        function f_plotEpsitime(obj)
            figure
            ax(1)=axes;

            plot(ax(1),obj.epsi.time_s(1:end-1),diff(obj.epsi.time_s),'k','LineWidth',obj.plot_properties.LineWidth)

            title(ax(1),sprintf('%s-%s-%s Epsi time',obj.Meta_Data.mission,obj.Meta_Data.vehicle_name,obj.Meta_Data.deployment))
            ylabel(ax(1),'\Delta time [sec]','Fontsize',obj.plot_properties.FontSize,'FontName',obj.plot_properties.FontName)
            xlabel(ax(1),'time [sec]','Fontsize',obj.plot_properties.FontSize,'FontName',obj.plot_properties.FontName)
            set(ax(1),'Fontsize',obj.plot_properties.FontSize,'FontName',obj.plot_properties.FontName)
            linkaxes(ax,'x')
        end
        %         function f_createProfiles(obj)
        %             EPSI_create_profiles_v2(obj.Meta_Data,...
        %                 obj.Meta_Data.PROCESS.Prmin_prof,...
        %                 obj.Meta_Data.PROCESS.Prcrit_prof,...
        %                 obj.Meta_Data.PROCESS.userConfirmsProfiles)
        %
        %             load(fullfile(obj.Meta_Data.paths.profiles,['Profiles_' obj.Meta_Data.deployment]));
        %
        %             switch obj.Meta_Data.vehicle_name
        %                 case 'FISH'
        %                     datachoice = 'datadown';
        %                     idxchoice = 'down';
        %                 case 'WW'
        %                     datachoice = 'dataup';
        %                     idxchoice = 'up';
        %                 otherwise
        %                     datachoice = 'datadown';
        %                     idxchoice = 'down';
        %             end
        %
        %             for iProf=1:length(CTDProfiles.(datachoice))
        %                 EPSI = EpsiProfiles.(datachoice){iProf};
        %                 CTD = CTDProfiles.(datachoice){iProf};
        %                 Profile = mod_epsilometer_merge_profile(obj.Meta_Data,CTD,EPSI);
        %                 Profile.profNum = iProf;
        %                 save(fullfile(obj.Meta_Data.paths.profiles,sprintf('Profile%03.0f',iProf)),'Profile','-v7.3');
        %             end
        %         end

        %         function obj = f_getProfileIndices(obj)
        %             % USAGE
        %             %   obj.f_getProfileIndices
        %             %
        %             % - Epsi_MakeMatFromRaw saves a merged pressure timeseries along
        %             % with individual .mat files and Epsi_MATFile_TimeIndex.
        %             % - This step takes the merged Epsi_PressureTimeseries and
        %             % separates into downcasts and upcasts (still to do)
        %             % - Saves Epsi_PressureTimeseries.mat in Meta_Data.paths.mat_data
        %             MATpath = obj.Meta_Data.paths.mat_data;
        %             CTD = EPSI_get_CTD_profile_indices(MATpath);
        %             obj = CTD;
        %             fprintf('... Epsi_PressureTimeseries saved in %s \n',MATpath)
        %         end

        function obj = f_calibrateTemperature(obj)
            % USAGE
            %   obj.Meta_Data = f_calibrateTemperature(obj);
            %
            if ~isfield(obj.Meta_Data.PROCESS,'nfft')
                error('Meta_Data.PROCESS.nfft is not defined. Add Meta_Data_Process info')
            end
            %             Meta_Data = obj.Meta_Data;
            %             save(fullfile(Meta_Data.paths.data,'Meta_Data'),'Meta_Data');

            %             try
            %                 load(fullfile(obj.Meta_Data.paths.profiles,['Profiles_' Meta_Data.deployment]));
            %             catch err
            %                 error('You need to create profiles before calibrating temperature for this deployment')
            %             end
            %
            switch obj.Meta_Data.vehicle_name
                case 'FISH'
                    obj.Meta_Data.PROCESS.profile_dir = 'down';
                    datachoice = 'datadown';
                    idxchoice = 'down';
                case {'WW','SEACYCLER'}
                    obj.Meta_Data.PROCESS.profile_dir = 'up';
                    datachoice = 'dataup';
                    idxchoice = 'up';
                otherwise
                    obj.Meta_Data.PROCESS.profile_dir = 'down';
                    datachoice = 'datadown';
                    idxchoice = 'down';
            end


            % NC 10/12/21 - Instead of using the longest profile, wait
            % until after profiles are created and define dTdV according to
            % calibrate_dTdV.m
            dTdV_process = 'old';
            switch dTdV_process
                case 'new'
                    obj.Meta_Data = process_calibrate_dTdV(obj.Meta_Data);
                case 'old'
                    % Load the pressure timeseries and find the downcast or upcast with the
                    % greatest range in pressure.
                    load(fullfile(obj.Meta_Data.paths.mat_data,'PressureTimeseries.mat'));
                    %             switch datachoice
                    %                 case 'dataup'
                    %                     profLengths = PressureTimeseries.endup-PressureTimeseries.startup;
                    %                     pRange = PressureTimeseries.P(PressureTimeseries.endup) -...
                    %                                 PressureTimeseries.P(PressureTimeseries.startup);
                    %                 case 'datadown'
                    %                     profLengths = PressureTimeseries.enddown-PressureTimeseries.startdown;
                    %                     pRange = PressureTimeseries.P(PressureTimeseries.enddown) -...
                    %                                 PressureTimeseries.P(PressureTimeseries.startdown);
                    %             end
                    profLengths = PressureTimeseries.endprof-PressureTimeseries.startprof;
                    pRange = PressureTimeseries.P(PressureTimeseries.endprof) -...
                        PressureTimeseries.P(PressureTimeseries.startprof);
                    % Find the longest profile
                    [~,idxProf] = max(pRange);

                    % Get (and merge if necessary) .mat data for this profile
                    tMin = PressureTimeseries.dnum(PressureTimeseries.startprof(idxProf));
                    tMax = PressureTimeseries.dnum(PressureTimeseries.endprof(idxProf));
                    Profile = obj.f_cropTimeseries(tMin,tMax);

                    %                 load(fullfile(obj.Meta_Data.paths.profiles,sprintf('Profile%03.0f',idxProf)));
                    %                 try
                    %                     Fs=obj.Meta_Data.AFE.FS;
                    %                 catch
                    %                     Fs=obj.Meta_Data.PROCESS.Fs; %Check this,  I think it's wrong
                    %                 end
                    %                 tscanLongestProfile = floor(0.8*length(Profile.epsi.time_s)/Fs);
                    %                 tscanDefault = 50;
                    %                 tscan = min([tscanDefault,tscanLongestProfile]);

                    % This is what actually calculates dTdV!
                    obj.Meta_Data=mod_epsi_temperature_spectra_v4(obj.Meta_Data,Profile);
            end

        end
        function obj = f_computeTurbulence(obj,Profile_or_profNum,saveData)
            % obj = f_computeTurbulence(obj,Profile_or_profNum,saveData)
            % Compute turbulence quantities (epsilon, chi, etc)
            %
            % USAGE
            %   Profile = f_computeTurbulence(obj,Profile,saveData);
            %
            % INPUTS
            %   Profile_or_profNum = Profile structure or number of profile
            %   already saved as a structure in L1
            %
            % OUTPUT
            %   Profile = structure similar to input Profile structure, but
            %   now with turbulence quantities added
            process_all_profiles = 0;
            if nargin<3
                saveData = 1;
                if nargin<2
                    process_all_profiles = 1;
                end
            end
            if ~any([isfield(obj.Meta_Data.PROCESS,'nfft'),isclassfield(obj.Meta_Data.PROCESS,'nfft')])
                obj.Meta_Data= obj.f_read_MetaProcess();
            end
            Meta_Data = obj.Meta_Data;
            if ~process_all_profiles
                obj = mod_epsilometer_calc_turbulence_v2(Meta_Data,Profile_or_profNum,saveData);
            elseif process_all_profiles
                profile_list = dir(fullfile(Meta_Data.paths.profiles,'Profile*.mat'));
                for p=1:length(profile_list)
                    fprintf('Building Profile%03.0f of %03.0f\n',p,length(profile_list))
                    obj = mod_epsilometer_calc_turbulence_v2(Meta_Data,p,saveData);
                end
            end
        end
        function obj = f_processNewProfiles(obj,varargin)
            % obj = f_processNewProfiles(obj,varargin)
            %
            % OPTIONAL ARGUMENTS
            %   'grid',P - 'grid' flag to grid some of the profile variables onto a
            %               standard pressure grid
            %            - P = the pressure array to use

            % If there are not yet any temperature calibration values,
            % calibrate the temperature probes.
            if obj.Meta_Data.AFE.t1.cal==0 && obj.Meta_Data.AFE.t2.cal==0
                obj = f_calibrateTemperature(obj);
            end
            if nargin>1
                obj = epsiProcess_processNewProfiles(obj,varargin);
            else
                obj = epsiProcess_processNewProfiles(obj);
            end
        end
        function obj = f_makeNewProfiles(obj)
            % obj = f_makeNewProfiles(obj)
            %
            % Makes new profiles but does not compute turbulence variables
            obj = epsiProcess_makeNewProfiles(obj);
        end
        function obj = f_interpolateProfileToP(obj,Profile,z)
            % obj = f_interpolateProfileToP(Profile,P)
            %
            % INPUTS:
            %   Profile - Profile structure with turbulence variables
            %   P - pressure array
            griddedProfile = epsiProcess_interpolate_Profile_to_P(Profile,z);
            obj = griddedProfile;
        end
        function obj = f_gridProfiles(obj,z)
            obj = epsiProcess_gridProfiles(obj,z);
        end
        function obj = f_cropTimeseries(obj,tMin,tMax)
            % Get a piece of timeseries structure that you can use to compute
            % turbulence variables.
            %
            % USAGE:
            %   Timeseries = ec.f_cropTimeseries(tMin,tMax)
            %
            % INPUTS:
            %   obj - epsi_class object
            %   tMin - minimum epsi.time_s (seconds)
            %   tMax - maximum epsi.time_s (seconds)
            %
            % OUTPUT:
            %   Timeseries - structure of epsi and ctd data to process turbulence
            %   variables
            tRange = [tMin,tMax];
            Meta_Data = obj.Meta_Data;
            Timeseries = epsiProcess_crop_timeseries(Meta_Data,tRange);
            obj = Timeseries;
        end
        function f_plotProfileSpectra(obj,Profile,depth,saveFig)
            if nargin<3
                saveFig=0;
            end
            plot_profile_and_spectra(Profile,depth,saveFig)
        end
        function f_clearProcessedData(obj)
            delete(fullfile(obj.Meta_Data.paths.mat_data,'*.mat'))
            delete(fullfile(obj.Meta_Data.paths.profiles,'*.mat'))
            delete(fullfile(obj.Meta_Data.paths.data,'Meta_Data.mat'));
        end
        function f_check_spectra(obj,id_profile)
            if nargin<2
                warning('f_check_spectra(id_profile)');
            end
            check_spectra(id_profile);
        end
        function f_check_grid_spectra(obj)
            check_grid_spectra();
        end
    end %end methods
end %end classdef
