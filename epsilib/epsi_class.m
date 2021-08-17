classdef epsi_class < handle
    %
    % epsi_class is a class of epsi data, processing functions, and
    % plotting functions
    %
    % Usage:
    % To create an epsi_class structure, you must first be within a
    % directory containing 1) an epsi config file, and 2) a directory
    % called 'raw' containing raw .ascii files
    %   >> obj = epsi_class(lastOrAll)
    %
    % INPUTS:i
    %   lastOrAll: 1=(default) loads last (most recent) file,
    %              2=loads all data in deployment
    
    properties
        Meta_Data %read from config file
        filename
        epsi
        ctd
        alt
        plot_properties
    end
    methods
        function obj=epsi_class(lastOrAll)
            if nargin<1
                lastOrAll=1;
            end
            %           TODO only define Meta data if epsi.mat and ctd.mat does not
            %           exist
            
            % Check to see if Meta_Data is already defined
            checkMD = exist('Meta_Data.mat','file');
            
            currDir = pwd;
            
            repeat = 1; %NC Initialize repeat flag to use if Meta_Data path names were not made on this machine
            while repeat==1
                switch checkMD
                    case 2 %Meta_Data, epsi, and ctd .mat files already exist
                        
                        fprintf('Initializing epsi_class with previously created Meta_Data \n')
                        load('Meta_Data.mat')
                        obj.Meta_Data = Meta_Data;
                        
                        % NC - If Meta_Data exists, check to see if path to mat data is
                        % available. If not, epsi_class was probably last generated on a
                        % different machine. Create a new Meta_Data.
                        checkMat = exist(fullfile(obj.Meta_Data.MATpath,'Epsi_MATfile_TimeIndex.mat'),'file');
                        if checkMat==0
                            checkMD = 0;
                        else
                            repeat = 0;
                        end
                        
                    case 0 %Meta_Data, epsi, and ctd .mat files do not already exist
                        
                        fprintf('Initializing epsi_class and creating new Meta_Data \n')
                        setupfile=dir('*config*');
                        repeat = 0;
                        
                        try
                            setup=mod_som_read_setup(setupfile.name);
                        catch
                            try
                                setupfile=dir('*_raw*');
                                setup=mod_som_read_setup_from_raw(setupfile(1).name);
                            catch
                                disp('mod_som_read_setup failed')
                            end
                        end
                        
                        try
                            obj.Meta_Data = fill_meta_data(setup);
                        catch
                            disp('fill_meta_data failed')
                            cd(currDir)
                        end
                end
            end
            
            % NC - check for s1 and s2 cal values. If they're
            % 0, manually input all probe numbers
            if obj.Meta_Data.AFE.s1.cal==0 || obj.Meta_Data.AFE.s2.cal==0
                obj.Meta_Data = obj.f_getSNshear;
                obj.Meta_Data = obj.f_getSNtemp;
            end
            
            % Read PROCESS Meta_Data from default text file
            obj.f_read_MetaProcess;
            
            % Define filename as path to raw data
            obj.filename=obj.Meta_Data.RAWpath;
            
            % Define plot properties
            obj.plot_properties = obj.f_getPlotProperties;
            
            % NC - Always check for new data by calling f_readData and then
            % load data into epsi class
            obj.f_readData();
            if lastOrAll==1
                obj = obj.f_getLastData();
            elseif lastOrAll==2
                obj = obj.f_getAllData();
            end
            cd(obj.Meta_Data.datapath)
        end
        function obj=f_readData(obj)
            dirs = {obj.Meta_Data.RAWpath; obj.Meta_Data.RAWpath; obj.Meta_Data.MATpath};
            % Convert ascii to mat
            Epsi_MakeMatFromRaw(dirs,obj.Meta_Data);
        end
        function obj = f_mergeData(obj)
            [epsi,ctd,alt] = merge_new_mat_files(obj.Meta_Data);
        end
        function obj = f_getLastData(obj);
            obj.f_readData;
            obj.epsi = obj.f_getLastEpsi();
            obj.ctd = obj.f_getLastCtd();
            obj.alt = obj.f_getLastAlt();
        end
        function obj = f_getAllData(obj);
            obj.f_mergeData;
            obj.epsi = obj.f_getAllEpsi();
            obj.ctd = obj.f_getAllCtd();
            obj.alt = obj.f_getAllAlt();
        end
        function obj=f_getLastEpsi(obj)
            load(fullfile(obj.Meta_Data.MATpath,'Epsi_MATfile_TimeIndex'));
            [~,idxLast] = max(Epsi_MATfile_TimeIndex.timeEnd);
            data = load(fullfile(obj.Meta_Data.MATpath,Epsi_MATfile_TimeIndex.filenames{idxLast}),'epsi');
            if isstruct(data.epsi)
                obj=data.epsi;
            else
                obj=[];
            end
        end
        function obj=f_getLastCtd(obj)
            load(fullfile(obj.Meta_Data.MATpath,'Epsi_MATfile_TimeIndex'));
            [~,idxLast] = max(Epsi_MATfile_TimeIndex.timeEnd);
            data = load(fullfile(obj.Meta_Data.MATpath,Epsi_MATfile_TimeIndex.filenames{idxLast}),'ctd');
            if isstruct(data.ctd)
                obj=data.ctd;
            else
                obj=[];
            end
        end
        function obj=f_getLastAlt(obj)
            load(fullfile(obj.Meta_Data.MATpath,'Epsi_MATfile_TimeIndex'));
            [~,idxLast] = max(Epsi_MATfile_TimeIndex.timeEnd);
            data = load(fullfile(obj.Meta_Data.MATpath,Epsi_MATfile_TimeIndex.filenames{idxLast}),'alt');
            if isstruct(data.alt)
                obj=data.alt;
            else
                obj=[];
            end
        end
        function obj=f_getAllEpsi(obj)
            data=load(fullfile(obj.Meta_Data.Epsipath,['epsi_' obj.Meta_Data.deployment '.mat']));
            if isstruct(data.epsi)
                obj=data.epsi;
            else
                obj=[];
            end
        end
        function obj=f_getAllCtd(obj)
            %obj.f_mergeData;
            if exist(fullfile(obj.Meta_Data.CTDpath,['ctd_' obj.Meta_Data.deployment '.mat'])) == 2
                data=load(fullfile(obj.Meta_Data.CTDpath,['ctd_' obj.Meta_Data.deployment '.mat']));
                if isstruct(data.ctd)
                    obj=data.ctd;
                else
                    obj=[];
                end
            else
                obj = [];
            end
        end
        function obj=f_getAllAlt(obj)
            %obj.f_mergeData;
            if exist(fullfile(obj.Meta_Data.CTDpath,['alt_' obj.Meta_Data.deployment '.mat'])) == 2
                data=load(fullfile(obj.Meta_Data.CTDpath,['alt_' obj.Meta_Data.deployment '.mat']));
                if isstruct(data.alt)
                    obj=data.alt;
                else
                    obj=[];
                end
            else
                obj = [];
            end
        end
        function obj=f_getPlotProperties(obj)
            % Set default plot properties (fonts, colors, sizes)
            obj = set_epsi_plot_properties;
        end
        function obj=f_getSNshear(obj)
            % Input shear probe serial numbers into Meta_Data and get
            % latest calibration values (Sv)
            %
            % USAGE
            %   obj.Meta_Data = f_getSNshear(obj)
            
            %             Meta_Data  = obj.Meta_Data;
            obj = set_SN_shear(obj.Meta_Data);
        end
        function obj=f_getSNtemp(obj)
            % Input temperature probe serial numbers into Meta_Data and get
            % latest calibration values (Sv)
            %
            % USAGE
            %   obj.Meta_Data = f_getSNtemp(obj)
            
            %             Meta_Data  = obj.Meta_Data;
            obj = set_SN_temp(obj.Meta_Data);
        end
        function obj=f_checkEpsiTime(obj,saveFig)
            if nargin<2
                saveFig=0;
            end
            check_epsi_time(obj,saveFig)
        end
        function f_plotCtd(obj)
            figure
            ax(1)=subplot(311);
            ax(2)=subplot(312);
            ax(3)=subplot(313);
            
            plot(ax(1),obj.ctd.ctdtime,obj.ctd.P,'k','LineWidth',obj.plot_properties.LineWidth)
            plot(ax(2),obj.ctd.ctdtime,obj.ctd.T,'r','LineWidth',obj.plot_properties.LineWidth)
            plot(ax(3),obj.ctd.ctdtime,obj.ctd.S,'b','LineWidth',obj.plot_properties.LineWidth)
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
        function f_plotEpsi(obj)
            figure
            ax(1)=subplot(311);
            ax(2)=subplot(312);
            ax(3)=subplot(313);
            
            plot(ax(1),obj.epsi.epsitime,detrend(obj.epsi.t1_volt,'constant'),'k','LineWidth',obj.plot_properties.LineWidth)
            hold(ax(1),'on')
            plot(ax(1),obj.epsi.epsitime,detrend(obj.epsi.t2_volt,'constant'),'m','LineWidth',obj.plot_properties.LineWidth)
            hold(ax(1),'on')
            
            plot(ax(2),obj.epsi.epsitime,detrend(obj.epsi.s1_volt,'constant'),'k','LineWidth',obj.plot_properties.LineWidth)
            hold(ax(2),'on')
            plot(ax(2),obj.epsi.epsitime,detrend(obj.epsi.s2_volt,'constant'),'m','LineWidth',obj.plot_properties.LineWidth)
            hold(ax(2),'on')
            
            plot(ax(3),obj.epsi.epsitime,detrend(obj.epsi.a1_g,'constant'),'k','LineWidth',obj.plot_properties.LineWidth)
            hold(ax(3),'on')
            plot(ax(3),obj.epsi.epsitime,detrend(obj.epsi.a2_g,'constant'),'m','LineWidth',obj.plot_properties.LineWidth)
            plot(ax(3),obj.epsi.epsitime,detrend(obj.epsi.a3_g,'constant'),'c','LineWidth',obj.plot_properties.LineWidth)
            hold(ax(3),'on')
            %ALB mean in legend
            legend(ax(1),{sprintf('t1 %3.2f V',nanmean(obj.epsi.t1_volt)),...
                sprintf('t2 %3.2f V',nanmean(obj.epsi.t2_volt))})
            legend(ax(2),{sprintf('s1 %3.2f V',nanmean(obj.epsi.s1_volt)),...
                sprintf('s2 %3.2f V',nanmean(obj.epsi.s2_volt))})
            legend(ax(3),{sprintf('a1 %3.2f g',nanmean(obj.epsi.a1_g)),...
                sprintf('a2 %3.2f g',nanmean(obj.epsi.a2_g)),...
                sprintf('a3 %3.2f g',nanmean(obj.epsi.a3_g))})
            
            ax(1).XTickLabel='';
            ax(2).XTickLabel='';
            title(ax(1),sprintf('%s-%s-%s',strrep(obj.Meta_Data.mission,'_','\_'),...
                strrep(obj.Meta_Data.vehicle_name,'_','\_'),...
                strrep(obj.Meta_Data.deployment,'_','\_')));
            ylabel(ax(1),'FPO7 [Volt]','Fontsize',obj.plot_properties.FontSize,'FontName',obj.plot_properties.FontName)
            ylabel(ax(2),'Shear [Volt]','Fontsize',obj.plot_properties.FontSize,'FontName',obj.plot_properties.FontName)
            ylabel(ax(3),'Accel [g]','Fontsize',obj.plot_properties.FontSize,'FontName',obj.plot_properties.FontName)
            ax(3).XLabel.String = 'seconds';
            set(ax(1),'Fontsize',obj.plot_properties.FontSize,'FontName',obj.plot_properties.FontName)
            set(ax(2),'Fontsize',obj.plot_properties.FontSize,'FontName',obj.plot_properties.FontName)
            set(ax(3),'Fontsize',obj.plot_properties.FontSize,'FontName',obj.plot_properties.FontName)
            linkaxes(ax,'x')
            
        end
        function f_plotAlti(obj)
            figure
            plotyy(obj.alt.alttime,obj.alt.dst,obj.ctd.ctdtime,obj.ctd.P);
        end
        function f_plotFallSpeed(obj)
            figure
            plot(obj.ctd.ctdtime,obj.ctd.dPdt,'.')
            hold on
            plot(obj.ctd.ctdtime,movmean(obj.ctd.dPdt,20),'.')
            xlabel('ctdtime (seconds)')
            ylabel('dPdt')
        end
        function f_plotEpsiCtd(obj,saveFig)
            if nargin<2
                saveFig=0;
            end
            plot_epsi_ctd(obj,saveFig)
        end
        function  [P11,f,noise]=f_calibrateEpsi(obj,tmid,tscan,makeFig,saveFig)
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
            if nargin<4
                makeFig=1;
                saveFig=0;
            end
            if nargin==4
                saveFig = 1;
            end
            [P11,f,noise] = mod_som_calibrate_epsi_tMid(obj,tmid,tscan,makeFig,saveFig);
        end
        function f_plotEpsitime(obj)
            figure
            ax(1)=axes;
            
            plot(ax(1),obj.epsi.epsitime(1:end-1),diff(obj.epsi.epsitime),'k','LineWidth',obj.plot_properties.LineWidth)
            
            title(ax(1),sprintf('%s-%s-%s Epsi time',obj.Meta_Data.mission,obj.Meta_Data.vehicle_name,obj.Meta_Data.deployment))
            ylabel(ax(1),'\Delta time [sec]','Fontsize',obj.plot_properties.FontSize,'FontName',obj.plot_properties.FontName)
            xlabel(ax(1),'time [sec]','Fontsize',obj.plot_properties.FontSize,'FontName',obj.plot_properties.FontName)
            set(ax(1),'Fontsize',obj.plot_properties.FontSize,'FontName',obj.plot_properties.FontName)
            linkaxes(ax,'x')
        end
        function f_createProfiles(obj)
            EPSI_create_profiles_v2(obj.Meta_Data,...
                obj.Meta_Data.PROCESS.Prmin_prof,...
                obj.Meta_Data.PROCESS.Prcrit_prof,...
                obj.Meta_Data.PROCESS.userConfirmsProfiles)
            
            load(fullfile(obj.Meta_Data.L1path,['Profiles_' obj.Meta_Data.deployment]));
            
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
            
            for iProf=1:length(CTDProfiles.(datachoice))
                EPSI = EpsiProfiles.(datachoice){iProf};
                CTD = CTDProfiles.(datachoice){iProf};
                Profile = mod_epsilometer_merge_profile(obj.Meta_Data,CTD,EPSI);
                Profile.profNum = iProf;
                save(fullfile(obj.Meta_Data.L1path,sprintf('Profile%03.0f',iProf)),'Profile','-v7.3');
            end
        end
        function obj = f_calibrateTemperature(obj,Profile)
            % USAGE
            %   obj.Meta_Data = f_calibrateTemperature(obj,Profile);
            %
            if ~isfield(obj.Meta_Data.PROCESS,'nfft')
                obj.Meta_Data= obj.f_read_MetaProcess();
            end
            Meta_Data = obj.Meta_Data;
            
            try
                load(fullfile(obj.Meta_Data.L1path,['Profiles_' Meta_Data.deployment]));
            catch err
                error('You need to create profiles before calibrating temperature for this deployment')
            end
            
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
            
            % For the temperature correction, we want a good chunk of data. Try
            % tscan = 50 or if the profile is too short, find the longest profile
            % in the deployment and set tscan to about 0.5-0.8 times the profile
            % length
            
            % Use the longest profile to do the temperature calibration
            [~,id_profile] = max(cellfun(@(C) length(C), EpsiProfiles.(idxchoice)));
            tscanAlternate = floor(0.8*length(EpsiProfiles.(datachoice){id_profile}.epsitime)/320);
            tscanDefault = 50;
            tscan = min([tscanDefault,tscanAlternate]);
            
            display=0;
            titleStr = strrep([obj.Meta_Data.mission ' ' obj.Meta_Data.vehicle_name ' ' obj.Meta_Data.deployment],'_','\_');
            
            obj.Meta_Data=mod_epsi_temperature_spectra_v2(obj.Meta_Data, ...
                EpsiProfiles.(datachoice){id_profile}, ...
                CTDProfiles.(datachoice){id_profile},...
                titleStr,id_profile,display,tscan);
                        
        end
        function obj = f_computeTurbulence(obj,Profile_or_profNum,saveData)
            % Compute turbulence quantities (epsilon, chi, etc)
            %
            % USAGE
            %   Profile = f_computeTurbulence(obj,Profile);
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
            if ~isfield(obj.Meta_Data.PROCESS,'nfft')
                obj.Meta_Data= obj.f_read_MetaProcess();
            end
            Meta_Data = obj.Meta_Data;
            obj = mod_epsilometer_calc_turbulence_v2(Meta_Data,Profile_or_profNum,saveData);
        end
        function obj = f_cropTimeseries(obj,tMin,tMax)
            % Get a short piece of timeseries structure that you can use to compute
            % turbulence variables.
            %
            % USAGE:
            %   Timeseries = ec.f_cropTimeseries(tMin,tMax)
            %
            % INPUTS:
            %   obj - epsi_class object
            %   tMin - minimum epsitime (seconds)
            %   tMax - maximum epsitime (seconds)
            %
            % OUTPUT:
            %   Timeseries - structure of epsi and ctd data to process turbulence
            %   variables
            tRange = [tMin,tMax];
            Meta_Data = obj.Meta_Data;
            Timeseries = crop_timeseries(Meta_Data,tRange);
            obj = Timeseries;
        end
        function f_plotProfileSpectra(obj,Profile,depth,saveFig)
            if nargin<2
                saveFig=0;
            end
            plot_profile_and_spectra(Profile,depth,saveFig)
        end
        function f_clearRawData(obj)
            delete(fullfile(obj.Meta_Data.CTDpath,'*.mat'))
            delete(fullfile(obj.Meta_Data.CTDpath,'*.mat'))
            delete(fullfile(obj.Meta_Data.Epsipath,'*.mat'))
            delete(fullfile(obj.Meta_Data.MATpath,'*.mat'))
            delete(fullfile(obj.Meta_Data.datapath,'Meta_Data.mat'));
        end
        function f_plotShadeFiles(obj,timeOrSamplenum)
            if nargin==1
                fprintf('Choose time or samplenum for x-axis\n')
            end
            epsi = obj.epsi;
            shade_different_files(obj,timeOrSamplenum);
        end
        function obj=f_read_MetaProcess(obj,filename)
            if nargin==1
                filename=fullfile(obj.Meta_Data.processpath,'EPSILON',...
                    'Meta_Data_Process.txt');
            end
            obj = read_MetaProcess(obj.Meta_Data,filename);
        end
        
        function obj=f_realtime_spectra(obj)
            x=2;
            while x<10
                obj = obj.f_getLastData;
                obj.f_calibrateEpsi(obj.epsi.epsitime(end)-2,4);
                pause(2)
            end
        end
        function f_realtime_epsilon(obj)
            
        end
    end
end


