classdef epsi_class
    %
    % epsi_class is a class of epsi data, processing functions, and
    % plotting functions
    %
    % Usage:
    % To create an epsi_class structure, you must first be within a
    % directory containing 1) an epsi config file, and 2) a directory
    % called 'raw' containing raw .ascii files
    %   >> obj = epsi_class

    properties
        Meta_Data %read from config file
        filename
        epsi
        ctd
        alt
        plot_properties
    end
    methods
        function obj=epsi_class()
            %           TODO only define Meta data if epsi.mat and ctd.mat does not
            %           exist
            
            % Check to see if Meta_Data is already defined
            checkMD = exist('Meta_Data.mat','file');
            currDir = pwd;
            switch checkMD
                case 2 %Meta_Data, epsi, and ctd .mat files already exist
                    
                    fprintf('Initializing epsi_class with previously created Meta_Data \n')
                    load('Meta_Data.mat')
                    obj.Meta_Data = Meta_Data;
   
                    checkEpsi = exist(fullfile(obj.Meta_Data.Epsipath,['epsi_' obj.Meta_Data.deployment '.mat']),'file');
                    
                case 0 %Meta_Data, epsi, and ctd .mat files do not already exist
                    
                    fprintf('Initializing epsi_class and creating new Meta_Data \n')
                    setupfile=dir('*config*');
                    
                    try
                        setup=mod_som_read_setup(setupfile.name);
                    catch
                        disp('mod_som_read_setup failed')
                    end
                    
                    try 
                        obj.Meta_Data = fill_meta_data(setup);
                    catch
                        disp('fill_meta_data failed')
                        cd(currDir)
                    end
            end
            
%             % Change directory to raw data path and 
%             cd(obj.Meta_Data.RAWpath)
%             %list=dir("*.ascii");
            obj.filename=obj.Meta_Data.RAWpath;
            % NC moved this to inside fill_meta_data
%             %ALB copy the raw file inside the raw folder. 
%             list_rawfile=dir("*.ascii");
%             for f=1:length(list_rawfile)
%                 copyfile(fullfile(list_rawfile(f).folder, ...
%                                   list_rawfile(f).name),  ...
%                                   "./raw/");
%             end
            obj.plot_properties.FontName='Sans';
            obj.plot_properties.FontSize=16;
            obj.plot_properties.LineWidth=2;
            obj.plot_properties.Colors.a1 = [129 27 112]./255;
            obj.plot_properties.Colors.a2 = [235 64 61]./255;
            obj.plot_properties.Colors.a3 = [245 199 118]./255;
            obj.plot_properties.Colors.s1 = [60 134 76]./255;
            obj.plot_properties.Colors.s2 = [173 215 136]./255;
            obj.plot_properties.Colors.t1 = [29 78 140]./255;
            obj.plot_properties.Colors.t2 = [78 173 173]./255;
            obj.plot_properties.Colors.P = [0.2353    0.2157    0.5098];
            obj.plot_properties.Colors.dPdt = [0.4851    0.4644    0.7140];
            obj.plot_properties.Colors.T = [185 38 26]./255;
            obj.plot_properties.Colors.S = [0 0 0];
            

            % NC moved this CTD cal section to inside fill_meta_data.m
%             % TODO fix the SBE name bug I do not want (1:3)
%             switch obj.Meta_Data.CTD.name
%                 case{'SBE49','SBE','S49'}
%                     try
%                         obj.Meta_Data.CTD.cal=get_CalSBE(obj.Meta_Data.CTD.CALfile(obj.Meta_Data.CTD.CALpath,obj.Meta_Data.CTD.SN));
%                         if(strcmp(obj.Meta_Data.CTD.SN,'0000'))
%                             obj.Meta_Data.CTD.SN=input('SBE49 SN(e.g. 0237):','s');
%                             obj.Meta_Data.CTD.cal=get_CalSBE(obj.Meta_Data.CTD.CALfile(obj.Meta_Data.CTD.CALpath,obj.Meta_Data.CTD.SN));
%                         end
%                     catch
%                         obj.Meta_Data.CTD.SN=input('SBE49 SN(e.g. 0237):','s');
%                         obj.Meta_Data.CTD.cal=get_CalSBE(obj.Meta_Data.CTD.CALfile(obj.Meta_Data.CTD.CALpath,obj.Meta_Data.CTD.SN));
%                     end
%                 case{'SBE41'}
%             end

            
            % Load epsi and ctd timeseries data into epsi class
                if isfile(fullfile(obj.Meta_Data.Epsipath,['epsi_' obj.Meta_Data.deployment '.mat']))
                    obj.epsi=obj.f_getEpsi();
                    obj.ctd=obj.f_getCtd();
                    obj.alt=obj.f_getAlt();
                else
                    obj.f_readData();
                    obj.epsi=obj.f_getEpsi();
                    obj.ctd=obj.f_getCtd();
                    obj.alt=obj.f_getAlt();
                end
            
            cd(obj.Meta_Data.datapath)
        end
        
        function obj=f_readData(obj)
            % mod_som_read_epsi_files(obj.filename,obj.Meta_Data);
            %addpath /Volumes/GoogleDrive/'Shared drives'/MOD-data-Epsilometer/Library/EPSILOMETER/EPSILON/toolbox/seawater2/
            mod_som_read_epsi_files_v2(obj.filename,obj.Meta_Data);
        end
        function obj=f_getEpsi(obj)
            data=load(fullfile(obj.Meta_Data.Epsipath,['epsi_' obj.Meta_Data.deployment '.mat']));
            if isstruct(data.epsi)
                obj=data.epsi;
             else
                obj=[];
            end
        end
        function obj=f_getCtd(obj)
            data=load(fullfile(obj.Meta_Data.CTDpath,['ctd_' obj.Meta_Data.deployment '.mat']));
            if isstruct(data.ctd)
                obj=data.ctd;
            else
                obj=[];
            end
        end
        function obj=f_getAlt(obj)
            data=load(fullfile(obj.Meta_Data.CTDpath,['alt_' obj.Meta_Data.deployment '.mat']));
            if isstruct(data.alt)
                obj=data.alt;
            else 
                obj=[];
            end
        end
        function obj=f_checkEpsiTime(obj)
            check_epsi_time(obj.Meta_Data)
        end
        function obj=f_checkCtdTime(obj)
            check_ctd_time(obj.Meta_Data)
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
        function f_plotEpsiCtd(obj)
            plot_epsi_ctd(obj)
        end
        function  [P11,f]=f_calibrateEpsi(obj,tmid,tscan,makeFig)
            % USAGE
            %   [P11,f] = obj.f_calibrateEpsi(tmid,tscan,makeFig)
            % 
            % INPUTS
            %   tmid = midpoint of scan (seconds)
            %   tscan = length of scan (seconds)
            %   makeFig = (OPTIONAL, flag to plot figure [0/1], default=1)
            %
            % OUTPUTS (OPTIONAL)
            %   P11 = structure of frequency spectra for each channel
            %   f = frequency array
            if nargin<4
                makeFig=1;
            end
            [P11,f] = mod_som_calibrate_epsi_tMid(obj,tmid,tscan,makeFig);
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
        end
        function f_computeTurbulence(obj,profileIdx)
            Meta_Data = obj.Meta_Data;
            if nargin==1
                 mod_epsilometer_batch_process_v2(Meta_Data)
            elseif nargin==2
            mod_epsilometer_batch_process_v2(Meta_Data,profileIdx)
            end
        end
        function f_clearRawData(obj)
            delete(fullfile(obj.Meta_Data.CTDpath,['ctd_' obj.Meta_Data.deployment '.mat']))
            delete(fullfile(obj.Meta_Data.CTDpath,['alt_' obj.Meta_Data.deployment '.mat']))
            delete(fullfile(obj.Meta_Data.Epsipath,['epsi_' obj.Meta_Data.deployment '.mat']))
            delete(fullfile(obj.Meta_Data.datapath,'Meta_Data.mat'));
        end   
        function f_plotShadeFiles(obj,timeOrSamplenum)
            if nargin==1
                fprintf('Choose time or samplenum for x-axis')
            end
            epsi = obj.epsi;
            shade_different_files(epsi,timeOrSamplenum);
        end
    end
end


