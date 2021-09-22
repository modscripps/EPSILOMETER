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
    % INPUTS:i
    %   lastOrAll: 1=(default) loads last (most recent) file,
    %              2=loads all data in deployment
    
    properties
        Meta_Data %read from config file
        filename
        plot_properties
        epsi
        ctd
        alt  
    end
    methods
        function obj=epsi_class(lastOrAll)
            if nargin<1
                lastOrAll=1;
            end
            %           TODO only define Meta data if epsi.mat and ctd.mat does not
            %           exist
            
            % Check to see if Meta_Data is already definedc
            checkMD = dir('Meta_Data.mat');
            
            repeat = 1; %NC Initialize repeat flag to use if Meta_Data path names were not made on this machine
            while repeat==1
                if ~isempty(checkMD) %Meta_Data already exists
                    
                    fprintf('Initializing epsi_class with previously created Meta_Data \n')
                    load('Meta_Data')
                    obj.Meta_Data = Meta_Data;
                    
                    % Check that processpath, datapath, and CALIpath are
                    % accessible. If not, Meta_Data was probably
                    % created on a different computer and you need to
                    % remake the paths.
                    % NC 9/22/21 - Always redefine the data path as the current
                    % directory
                    obj.Meta_Data.paths.data=pwd;
%                     idx_blank_space=strfind(obj.Meta_Data.paths.data,' ');
%                     
%                     if ~isempty(idx_blank_space)
%                         split_str=strsplit(Meta_Data.paths.data,' ');
%                         for i=1:length(split_str)-1
%                             split_str{i}=[split_str{i} '\ ']
%                         end
%                         Meta_Data.paths.data=[split_str{:}]
%                         obj.Meta_Data.paths.data=Meta_Data.paths.data;
%                     end
                    if  ~isdir(obj.Meta_Data.paths.process_library) || ...
                        ~isdir(obj.Meta_Data.paths.data) || ...
                        ~isdir(obj.Meta_Data.paths.calibration) || ...
                        ~isclassfield(obj.Meta_Data.paths,'raw_data')
%                         obj.Meta_Data = set_epsi_paths(obj.Meta_Data);
                        
                        % Find the epsi library and add it as process path
                        spltpath=strsplit(path,':');
                        epsilib_path=spltpath{~cellfun(@isempty, ...
                            cellfun(@(x) ...
                            strfind(x,'epsilib'),spltpath, ...
                            'UniformOutput',false))};
                        
                        obj.Meta_Data.paths.process_library=fileparts(epsilib_path);
                        obj.Meta_Data.paths.calibration = fullfile(obj.Meta_Data.paths.process_library,'CALIBRATION','ELECTRONICS');

                        % Set epsi paths and define suffix for raw files
                        obj.Meta_Data = epsiSetup_set_epsi_paths(obj.Meta_Data);
                        obj.Meta_Data = epsiSetup_get_raw_suffix(obj.Meta_Data);
                        
                        Meta_Data = obj.Meta_Data;
                        save(fullfile(obj.Meta_Data.paths.data,'Meta_Data'),'Meta_Data');
                    else
                        repeat = 0;
                        checkMD = [];
                    end
                    
                elseif isempty(checkMD) %Meta_Data, epsi, and ctd .mat files do not already exist
                    
                    fprintf('Initializing epsi_class and creating new Meta_Data \n')
                    repeat = 0;
                    
                    % Set epsi paths and define suffix for raw files
                    obj.Meta_Data = epsiSetup_set_epsi_paths(obj.Meta_Data);
                    obj.Meta_Data = epsiSetup_get_raw_suffix(obj.Meta_Data);
                    
                    % First, try reading configuration data from the
                    % file. If that doesn't work, try reading from a
                    % configuration file.
%                                         try
%                                             setupfile=dir(fullfile(obj.Meta_Data.paths.raw_data,['*' obj.Meta_Data.rawfileSuffix]));
% %                                             setup=mod_som_read_setup_from_raw(setupfile(1).name);
%                                         catch
%                 end
                    try
                        setupfile=dir('*config*');
                        setup=mod_som_read_setup_from_config(setupfile.name);
                    catch
                        error('mod_som_read_setup failed - do you have a config file?')
                    end
                    % end
                    
                    try
                        obj.Meta_Data = epsiSetup_fill_meta_data(setup);
                        
                        fprintf('Meta_Data.paths.process_library is %s \n',obj.Meta_Data.paths.process_library);
                        fprintf('Meta_Data.paths.data is %s \n',obj.Meta_Data.paths.data);
                    catch
                        error('fill_meta_data failed')
                    end
                    
                    % Set epsi paths and define suffix for raw files
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
            
            % NC - check for s1 and s2 cal values. If they're
            % 0, manually input all probe numbers
            if obj.Meta_Data.AFE.s1.cal==0 || obj.Meta_Data.AFE.s2.cal==0
                obj.Meta_Data = obj.f_getSNshear;
                obj.Meta_Data = obj.f_getSNtemp;
            end
            
            % Read PROCESS Meta_Data from default text file
            obj.f_read_MetaProcess;
            
            % Define filedir as path to raw data
            obj.filename=obj.Meta_Data.paths.raw_data;
            
            % Define plot properties
            obj.f_getPlotProperties;
            
            % Don't automatically read data. If you ran autoreadEpsi, you
            % probably don't need to and you 
%             %NC - Always check for new data by calling f_readData and then
%             %load data into epsi class
%             obj.f_readData();
%             if lastOrAll==1
%                 obj = obj.f_getLastData();
%             elseif lastOrAll==2
%                 obj = obj.f_getAllData();
%             end
%             cd(obj.Meta_Data.paths.data)
        end
        
        
        function obj=f_read_MetaProcess(obj,filename)
            if nargin==1
                filename=fullfile(obj.Meta_Data.paths.process_library,'Meta_Data_Process',...
                    'Meta_Data_Process.txt');
            end
            Meta_Data = epsiSetup_read_MetaProcess(obj.Meta_Data,filename);
            obj.Meta_Data = Meta_Data;
        end
        function obj=f_readData(obj,varargin)

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
            dirs = {obj.Meta_Data.paths.raw_data; obj.Meta_Data.paths.mat_data};
            epsiProcess_convert_new_raw_to_mat(dirs,obj.Meta_Data,'noSync');
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
            data = load(fullfile(obj.Meta_Data.paths.mat_data,TimeIndex.filenames{idxLast}),'ctd');
            
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
                data.ctd = epsiProcess_merge_mat_files(data.ctd,data2.ctd);
                data.alt = epsiProcess_merge_mat_files(data.alt,data2.alt);
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
        
        
        function  [P11,f,noise,ax]=f_plot_spectraAtTmid(obj,tmid,tscan,nSec,makeFig,saveFig,replaceData,ax)
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
            %   tmid = midpoint of scan (seconds)
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
                obj.Meta_Data= obj.f_read_MetaProcess();
            end
            Meta_Data = obj.Meta_Data;
            save(fullfile(Meta_Data.paths.data,'Meta_Data'),'Meta_Data');
            
%             try
%                 load(fullfile(obj.Meta_Data.paths.profiles,['Profiles_' Meta_Data.deployment]));
%             catch err
%                 error('You need to create profiles before calibrating temperature for this deployment')
%             end
%             
            switch obj.Meta_Data.vehicle_name
                case 'FISH'
                    datachoice = 'datadown';
                    idxchoice = 'down';
                case 'WW'
                    datachoice = 'dataup';
                    idxchoice = 'up';
                otherwise
                    datachoice = 'datadown';
                    idxchoice = 'down';
            end
            
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
                
                obj.Meta_Data=mod_epsi_temperature_spectra_v3(obj.Meta_Data,Profile);
            
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
            if nargin<3
                saveData = 1;
            end
            if ~any([isfield(obj.Meta_Data.PROCESS,'nfft'),isclassfield(obj.Meta_Data.PROCESS,'nfft')])
                obj.Meta_Data= obj.f_read_MetaProcess();
            end
            Meta_Data = obj.Meta_Data;
            obj = mod_epsilometer_calc_turbulence_v2(Meta_Data,Profile_or_profNum,saveData);
        end
        function obj = f_processNewProfiles(obj,varargin)
            % obj = f_processNewProfiles(obj,varargin)
            %
            % OPTIONAL ARGUMENTS
            %   'grid',P - 'grid' flag to grid some of the profile variables onto a
            %               standard pressure grid
            %            - P = the pressure array to use
            %
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
    end %end methods
end %end classdef


