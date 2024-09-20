function Meta_Data=epsiSetup_fill_meta_data(Meta_Data,setup)
%% fill up the Meta_Data structure.
% There are two places where Meta_Data come from:
%    setup = a structure that is read from the instrument or from a
%            configuration file
%    Meta_Data_process_file = a .txt file with process 
% This function only adds Meta_Data from the setup structure.
%
%  It is used to process and organize epsi data.
%  The user will all necessarry informations to understand the epsi data
%  - electronics name, serial numbers, revision
%  - epsi number of channel and their names
%  - information on the CTD
%  - information on the epsi vehicle (epsifish, wirewalker, glider,...)
%  - firmware revision
%
%  TODO look at the TODO in this file. Most of them imply changes in the firmware.
%
%  written by A. Le Boyer 10/22/2020

% get setup field name.
% used to fill up the Meta_Data fields
setup_fields=fieldnames(setup);

spltpath=strsplit(path,':');
% If epsilib directory is immediately inside EPSILOMETER directory, the
% next two lines should appropriately define the process directory
% ALB to NC: Pretty smart move :).
epsilib_path={spltpath{~cellfun(@isempty, ...
                               cellfun(@(x) ...
                               strfind(x,'epsilib'),spltpath, ...
                               'UniformOutput',false))}};
                           
% These paths should have already been defined before entering this
% function -DON'T USE PWD BECAUSE YOU MIGHT BE WORKING IN A DIFFERENT
% DIRECTORY THAN WHERE THE DATA ARE STORED
%Meta_Data.paths.process_library=fileparts(epsilib_path{cellfun(@length,epsilib_path)==min(cellfun(@length,epsilib_path))});
%Meta_Data.paths.data=pwd;

Meta_Data.mission=setup.mission_name;
Meta_Data.vehicle_name=setup.vehicle_name;
Meta_Data.deployment=setup.SDIO.prefix_file;

% Get start time
Meta_Data.start_dnum = setup.start_dnum;

% Add profile direction based on vehicle name %NC 3/1/22
switch lower(Meta_Data.vehicle_name)
    case {'fish'}
        Meta_Data.PROCESS.profile_dir = 'down';
    case {'ww','seacycler','apex'}
        Meta_Data.PROCESS.profile_dir = 'up';
    otherwise
        Meta_Data.PROCESS.profile_dir = 'down';
end

%% get controler (CTL) name
Controlernames={'SOM','MADRE','PERSISTOR'};% hard coded name of potential CONTROLER
wh_CTL=cellfun(@(y) find(cellfun(@(x) strcmp(x,y),Controlernames)),setup_fields,'un',0);
wh_CTL=Controlernames{wh_CTL{~cellfun(@isempty,wh_CTL)}};


Meta_Data.CTL.name=wh_CTL;
Meta_Data.CTL.rev='rev4'; %TODO get info from config file
Meta_Data.CTL.SN='000';      %TODO get info from config file

%% get analog front end (AFE) name
Analog_names={'EFE','MAP','FLUO'};% hard coded name of potential CONTROLER
% find the kind of CTD we used from the setup file.
wh_AFE=cellfun(@(y) find(cellfun(@(x) strcmp(x,y),Analog_names)),setup_fields,'un',0);
wh_AFE=Analog_names{wh_AFE{~cellfun(@isempty,wh_AFE)}};

Meta_Data.AFE.name=wh_AFE;
Meta_Data.AFE.rev='EFErev4'; %TODO get info from config file
Meta_Data.AFE.SN='000';      %TODO get info from config file

% Meta_Data.paths.calibration=fullfile('..','CALIBRATION',[Meta_Data.CTL.rev '_' ...
%     Meta_Data.CTL.SN '-' Meta_Data.AFE.rev '_' Meta_Data.AFE.SN]);
Meta_Data.paths.calibration=fullfile(Meta_Data.paths.process_library,'CALIBRATION','ELECTRONICS');
%% set process parameters
Meta_Data.PROCESS.nb_channels = setup.(wh_AFE).nb_channel;
Meta_Data.PROCESS.channels=cellfun(@(x) x.name, setup.(wh_AFE).sensors, 'un',0);
Meta_Data.PROCESS.recording_mod='SD';

% Make 'timeseries' from 'channels'. (Add _g to acceleration channels and

for n=1:numel(Meta_Data.PROCESS.channels)
    if contains(Meta_Data.PROCESS.channels{n},'a')
        Meta_Data.PROCESS.timeseries{n} = sprintf('%s_g',Meta_Data.PROCESS.channels{n});
    elseif contains(Meta_Data.PROCESS.channels(n),{'s','t'})
        Meta_Data.PROCESS.timeseries{n} = sprintf('%s_volt',Meta_Data.PROCESS.channels{n});
    elseif contains(Meta_Data.PROCESS.channels(n),'c')
        Meta_Data.PROCESS.timeseries{n} = sprintf('%s_count',Meta_Data.PROCESS.channels{n});
    end
end

Meta_Data.AFE.shearcal_path=fullfile(Meta_Data.paths.process_library,'CALIBRATION','SHEAR_PROBES');
Meta_Data.AFE.tempcal_path=fullfile(Meta_Data.paths.process_library,'CALIBRATION','FPO7');

for i=1:Meta_Data.PROCESS.nb_channels
    sensor=setup.(wh_AFE).sensors{i};
    wh_name=sensor.name;
    % for the temp circuit TDIFF is only used with Bipolar
    % but we are not using TDIFF
    Meta_Data.AFE.(wh_name).SN=sensor.sn.';
    Meta_Data.AFE.(wh_name).cal=sensor.cal;
    switch wh_name
        case {'t1','t2','s1','s2','f1','f2','c1','c2'}
            Meta_Data.AFE.(wh_name).full_range=2.5;
        case {'a1','a2','a3'}
            Meta_Data.AFE.(wh_name).full_range=1.8;
        otherwise
                Meta_Data.AFE.(wh_name).full_range=input("What is the full range of that channel");
    end

    switch sensor.register.CONFIG_0
        case '1E0'
            Meta_Data.AFE.(wh_name).ADCconf='Unipolar';
            Meta_Data.AFE.temp_circuit='none';
        case '9E0'
            Meta_Data.AFE.(wh_name).ADCconf='Bipolar';
            Meta_Data.AFE.temp_circuit='tdiff';
        otherwise
            switch wh_name
                case {'t1','t2','a1','a2','a3'}
                    Meta_Data.AFE.(wh_name).ADCconf='Unipolar';
                    Meta_Data.AFE.temp_circuit='none';
                case {'s1','s2'}
                    Meta_Data.AFE.(wh_name).ADCconf='Bipolar';
                    Meta_Data.AFE.temp_circuit='none';
                otherwise
                    disp("Wrong CONFIG_0");
            end

    end

    switch setup.EFE.sensors{i}.register.FILTER_0
        case '6003C'
            Meta_Data.AFE.(wh_name).ADCfilter='sinc4';
            Meta_Data.AFE.FS=320;
        otherwise
            disp('no FILTER_0 REGISTER set to default')
            Meta_Data.AFE.(wh_name).ADCfilter='sinc4';
            Meta_Data.AFE.FS=320;
    end

end
Meta_Data=mod_som_get_shear_probe_calibration_v2(Meta_Data);
<<<<<<< Updated upstream
=======
Meta_Data=mod_som_get_temp_probe_calibration(Meta_Data);

>>>>>>> Stashed changes
Meta_Data.AFE.shear='CAmp1.0'; %TODO get info from config file

Meta_Data.Firmware.version='mod_som_som_eferev3_sdio_sampling_app_07152020.sls'; %TODO get info from config file

%% add auxillary device field

%I am commenting the all CTD section It needs to be read from the .modraw


CTDnames={'SBE','SBE49','SBE41','SB41','RBR','S49','SB49'};% hard coded name of potential CTD we will use with epsi
setup_fields=fieldnames(setup);
% % find the kind of CTD we used from the setup file.
wh_CTD=cellfun(@(y) find(cellfun(@(x) strcmp(x,y),CTDnames)),setup_fields,'un',0);
wh_CTD=CTDnames{wh_CTD{~cellfun(@isempty,wh_CTD)}};
% 
% if isfield(setup,wh_CTD)
Meta_Data.CTD.name = setup.(wh_CTD).header.';
% %TODO add the serial number in the SBE49 setup file. Maybe I want to get that after a the ds cmd.
% % Also use SBE in TPS and NOT engineer format.
% Meta_Data.CTD.SN   = num2str(str2double(setup.(wh_CTD).sn),'%04.0f');
Meta_Data.CTD.sample_per_record   = setup.(wh_CTD).sample_data_per_record;
% Meta_Data.CTD.CALpath   = fullfile(Meta_Data.paths.process_library,'CALIBRATION','SBE49');
% 
% Meta_Data.CTD.CALfile   = @(x,y) fullfile(x,[y '.cal']);
% 
% % % for pressure tests
% % disp('fill_meta_data.m line 169 - !!! auto-fill ctd SN !!!')
% % Meta_Data.CTD.SN = '0537';
% % Meta_Data.CTD.cal = get_CalSBE(Meta_Data.CTD.CALfile(Meta_Data.CTD.CALpath,Meta_Data.CTD.SN));
% 
% % TO DO - Get CTD SN from raw file header
% % NC 13 Aug. 2022 - For now, SBE serial number is in the file header but
% % not in the SOM3 line. Open the setup file back up and read the serial
% % number from the header
% if isfield(setup,'filepath')
% 
%     fid=fopen(setup.filepath);
%     total_str = fread(fid,'*char');
%     [ind_sbe_start,ind_sbe_stop, ind_settings_tokens] = regexp(total_str.','\CTD.SerialNum([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end','tokenExtents');
% 
%     % A really crude way to get the serial number - might not work every
%     % time. For example, if there are spaces before the number. Brute force
%     % for now.
%     SBE_sn = total_str(ind_sbe_start+15:ind_sbe_start+18).';
% 
%     Meta_Data.CTD.name = 'SBE49';
%     Meta_Data.CTD.SN = SBE_sn;
% 
% end
% 
% switch Meta_Data.CTD.name
%     case{'SBE49','SBE','S49','SB49'}
%         try
%             Meta_Data.CTD.cal=get_CalSBE(Meta_Data.CTD.CALfile(Meta_Data.CTD.CALpath,Meta_Data.CTD.SN));
%             if(strcmp(Meta_Data.CTD.SN,'0000'))
%                 Meta_Data.CTD.SN=input('**** SBE49 SN (e.g. 0237):','s');
%                 pause(3)
%                 Meta_Data.CTD.cal=get_CalSBE(Meta_Data.CTD.CALfile(Meta_Data.CTD.CALpath,Meta_Data.CTD.SN));
%             end
%         catch
%             Meta_Data.CTD.SN=input('**** SBE49 SN (e.g. 0237):','s');
%             pause(3)
%             Meta_Data.CTD.cal=get_CalSBE(Meta_Data.CTD.CALfile(Meta_Data.CTD.CALpath,Meta_Data.CTD.SN));
%         end
%     case{'SBE41'}
% end
% end
% 
Meta_Data.SDIO=setup.SDIO;

% Meta_Data=mod_som_define_epsi_meta_data(Meta_Data);



% fprintf('Saving Meta_Data in datapath \n')
% save(fullfile(Meta_Data.paths.data,'Meta_Data.mat'),'Meta_Data');

% TODO: It does not work if the data does not have all the channels
% .     I ll change that if needed
