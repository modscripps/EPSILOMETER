function [data] = mod_som_read_epsi_files_v4(filename,Meta_Data)

% For now, only works for 'v4' of som_acq
% It only works with a single file, since we're always calling it from a
% loop  in Epsi_MakeMatFromRaw
%
% ALB: I added seg and spec structures:
% Seg is a structure with the segments compiled in
% the fill_segment task inside efeobp.c
% Spec is a structure with the spectra  compiled in
% the cpt_spectra task inside efeobp.c. The spectra should match the segment found in the seg structure.
%
% TO DO:
%   Add backwards compatability - look at mod_som_read_epsi_files_v3.m
%
% CALLED BY:
%   epsiProcess_convert_new_raw_to_mat
%
% INPUTS
%   filename
%   Meta_Data
%
% OUTPUTS
%   epsi  (intermediate processing in 'efe' structure)
%   ctd   (intermediate processing in 'sbe' structure)
%   alt   (intermediate processing in 'alti' structure)
%   isap   (intermediate processing in 'isap' structure)
%   act   (intermediate processing in 'actu' structure)
%   vnav  (intermediate processing in 'vecnav' structure)
%   gps   (intermediate processing in 'gpsmeta' structure
%   seg   (intermediate processing in 'seg' structure)
%   spec  (intermediate processing in 'spec' structure)
%   eco  (intermediate processing in 'eco' structure)
%   ttv  (intermediate processing in 'ttv' structure)
%
% 11/28/21 aleboyer@ucsd.edu  updated from v3
% Nicole Couto adapted from Arnaud LeBoyer's mod_som_read_epsi_raw.m
% June 2021
% 2024/08/08 aleboyer: adding isap struct
% -------------------------------------------------------------------------

% Define a constant for salinity calculation
c3515 = 42.914;

%% Open file and save contents as 'str'
fprintf("   Open %s \r\n",filename)
fid = fopen(filename);
fseek(fid,0,1);
frewind(fid);
str = fread(fid,'*char')';
fclose(fid);

% Get time now
Meta_Data.start_dnum = now;

%% Get indices and tokens for each data type you will process
% ind_*_start  = starting indices of all matches
% ind_*_end    = ending indices of all matches

[ind_som_start, ind_som_stop]            = regexp(str,'\$SOM3([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
[ind_dcal_start, ind_dcal_stop]          = regexp(str,'\$DCAL([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
[ind_gps_start      , ind_gps_stop]      = regexp(str,'\$GPGGA([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
if isempty(ind_gps_start)
    [ind_gps_start      , ind_gps_stop]  = regexp(str,'\$INGGA([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
end

[ind_efe_start, ind_efe_stop]            = regexp(str,'\$EFE([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
[ind_sbe_start, ind_sbe_stop]            = regexp(str,'\$SB49([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
if isempty(ind_sbe_start)
    [ind_sbe_start, ind_sbe_stop]        = regexp(str,'\$SB41([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
end
[ind_alt_start      , ind_alt_stop]      = regexp(str,'\$ALT([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
[ind_isap_start      , ind_isap_stop]      = regexp(str,'\$ISAP([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
[ind_act_start      , ind_act_stop]      = regexp(str,'\$ACTU([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
[ind_vnav_start     , ind_vnav_stop]     = regexp(str,'\$VNMAR([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
if isempty(ind_vnav_start)
    [ind_vnav_start     , ind_vnav_stop] = regexp(str,'\$VNYPR([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
end
[ind_seg_start      , ind_seg_stop]      = regexp(str,'\$SEGM([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
[ind_spec_start     , ind_spec_stop]     = regexp(str,'\$SPEC([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
[ind_avgspec_start  , ind_avgspec_stop]  = regexp(str,'\$AVGS([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
[ind_dissrate_start , ind_dissrate_stop] = regexp(str,'\$RATE([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
[ind_apf0_start     , ind_apf0_stop]     = regexp(str,'\$APF0([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
[ind_apf1_start     , ind_apf1_stop]     = regexp(str,'\$APF1([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
[ind_apf2_start     , ind_apf2_stop]     = regexp(str,'\$APF2([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
[ind_fluor_start    , ind_fluor_stop]    = regexp(str,'\$ECOP([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
[ind_ttv_start      , ind_ttv_stop]      = regexp(str,'\$TTV1([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');



%% Get experiment, cruise, vehicle, pressure case and fish flag from setup, 
newSetup_flag=contains(str,'CTD.experiment=');
if newSetup_flag
    experiment_str    = str(strfind(str,'CTD.experiment=')+(0:100));
    cruise_str        = str(strfind(str,'CTD.cruise=')+(0:100));
    vehicle_str       = str(strfind(str,'CTD.vehicle=')+(0:100));
    pressure_case_str = str(strfind(str,'CTD.fish_pc=')+(0:100));
    fishflag_str      = str(strfind(str,'CTD.fishflag=')+(0:100));

    experiment_str    = experiment_str(1:find(uint8(experiment_str)==10,1,'first'));
    cruise_str        = cruise_str(1:find(uint8(cruise_str)==10,1,'first'));
    vehicle_str       = vehicle_str(1:find(uint8(vehicle_str)==10,1,'first'));
    pressure_case_str = pressure_case_str(1:find(uint8(pressure_case_str)==10,1,'first'));
    fishflag_str      = fishflag_str(1:find(uint8(fishflag_str)==10,1,'first'));

    experiment_name    = strsplit(experiment_str,'=');
    cruise_name        = strsplit(cruise_str,'=');
    vehicle_name       = strsplit(vehicle_str,'=');
    pressure_case_name = strsplit(pressure_case_str,'=');
    fishflag_name      = strsplit(fishflag_str,'=');


    experiment_name    = experiment_name{2}(1:end-2);
    cruise_name        = cruise_name{2}(1:end-2);
    vehicle_name       = vehicle_name{2}(1:end-2);
    pressure_case_name = pressure_case_name{2}(1:end-2);
    fishflag_name      = fishflag_name{2}(1:end-2);
else
    experiment_name    = [];
    cruise_name        = [];
    vehicle_name       = [];
    pressure_case_name = [];
    fishflag_name      = [];

end

%% get CTD cal coef
header_length=strfind(str,'END_FCTD_HEADER_START_RUN');
str_SBEcalcoef_header=str(strfind(str,'SERIALNO'):header_length);
if ~isempty(str_SBEcalcoef_header)
    SBEcal=get_CalSBE_from_modraw_header(str_SBEcalcoef_header);
    Meta_Data.CTD.cal=SBEcal;
else
    if ~isempty(ind_sbe_start)
        fprintf("SBE cal coef are missing in file %s ",filename)
        listfile=dir(Meta_Data.paths.raw_data);
        list_fullfilename=fullfile({listfile.folder},{listfile.name});
        idx_file=find(cellfun(@(x) strcmp(x,filename),list_fullfilename));
        count=0;
        %ALB If SBEcal missing  ALB using the SBEcal from the previous 10 modraw
        %files
        while isempty(str_SBEcalcoef_header)
            count=count+1;
            fprintf("  Using Open %s \r\n",list_fullfilename{idx_file-count})
            fid1 = fopen(list_fullfilename{idx_file-count});
            fseek(fid1,0,1);
            frewind(fid1);
            str1 = fread(fid1,'*char')';
            fclose(fid1);
            % get CTD cal coef
            header_length=strfind(str1,'END_FCTD_HEADER_START_RUN');
            str_SBEcalcoef_header=str1(strfind(str1,'SERIALNO'):header_length);
            if ~isempty(str_SBEcalcoef_header)
                SBEcal=get_CalSBE_from_modraw_header(str_SBEcalcoef_header);
            end
            if count>10
                error("No SBE data in modraw %s. \n Check the file ",filename)
            end
        end
    end

end
%% get Epsi probe Serial Numbers from Header
str_EPSICHANNEL_start  = strfind(str,'Fish probes serial numbers');
str_EPSICHANNEL_end    = strfind(str,'end probes serial numbers');
str_EPSICHANNEL_header = str(str_EPSICHANNEL_start:str_EPSICHANNEL_end-8);
epsi_probes=[];
if ~isempty(str_EPSICHANNEL_header)
    epsi_probes=parse_epsi_channel_string(str_EPSICHANNEL_header);
    Meta_Data.AFE.t1=epsi_probes.ch1;
    Meta_Data.AFE.t2=epsi_probes.ch2;
    Meta_Data.AFE.s1=epsi_probes.ch3;
    Meta_Data.AFE.s2=epsi_probes.ch4;
else
    warning('No Epsi probe serial number in file %s',filename)

end


%% Get indices and tokens for each data type you will process
% ind_*_start  = starting indices of all matches
% ind_*_end    = ending indices of all matches

<<<<<<< Updated upstream
% [ind_som_start, ind_som_stop]            = regexp(str,'\$SOM3([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
% [ind_dcal_start, ind_dcal_stop]          = regexp(str,'\$DCAL([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
% [ind_gps_start      , ind_gps_stop]      = regexp(str,'\$GPGGA([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
% if isempty(ind_gps_start)
%     [ind_gps_start      , ind_gps_stop]  = regexp(str,'\$INGGA([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
% end
% 
% [ind_efe_start, ind_efe_stop]            = regexp(str,'\$EFE([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
% [ind_sbe_start, ind_sbe_stop]            = regexp(str,'\$SB49([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
% if isempty(ind_sbe_start)
%     [ind_sbe_start, ind_sbe_stop]        = regexp(str,'\$SB41([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
% end
% [ind_alt_start      , ind_alt_stop]      = regexp(str,'\$ALT([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
% [ind_isap_start      , ind_isap_stop]      = regexp(str,'\$ISAP([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
% [ind_act_start      , ind_act_stop]      = regexp(str,'\$ACTU([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
% [ind_vnav_start     , ind_vnav_stop]     = regexp(str,'\$VNMAR([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
% if isempty(ind_vnav_start)
%     [ind_vnav_start     , ind_vnav_stop] = regexp(str,'\$VNYPR([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
% end
% [ind_seg_start      , ind_seg_stop]      = regexp(str,'\$SEGM([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
% [ind_spec_start     , ind_spec_stop]     = regexp(str,'\$SPEC([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
% [ind_avgspec_start  , ind_avgspec_stop]  = regexp(str,'\$AVGS([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
% [ind_dissrate_start , ind_dissrate_stop] = regexp(str,'\$RATE([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
% [ind_apf0_start     , ind_apf0_stop]     = regexp(str,'\$APF0([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
% [ind_apf1_start     , ind_apf1_stop]     = regexp(str,'\$APF1([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
% [ind_apf2_start     , ind_apf2_stop]     = regexp(str,'\$APF2([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
% [ind_fluor_start    , ind_fluor_stop]    = regexp(str,'\$ECOP([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
% [ind_ttv_start      , ind_ttv_stop]      = regexp(str,'\$TTVP([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
=======
[ind_som_start, ind_som_stop]            = regexp(str,'\$SOM3([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
[ind_dcal_start, ind_dcal_stop]          = regexp(str,'\$DCAL([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
[ind_gps_start      , ind_gps_stop]      = regexp(str,'\$GPGGA([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
if isempty(ind_gps_start)
    [ind_gps_start      , ind_gps_stop]  = regexp(str,'\$INGGA([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
end

[ind_efe_start, ind_efe_stop]            = regexp(str,'\$EFE([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
[ind_sbe_start, ind_sbe_stop]            = regexp(str,'\$SB49([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
if isempty(ind_sbe_start)
    [ind_sbe_start, ind_sbe_stop]        = regexp(str,'\$SB41([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
end
[ind_alt_start      , ind_alt_stop]      = regexp(str,'\$ALT([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
[ind_isap_start      , ind_isap_stop]      = regexp(str,'\$ISAP([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
[ind_act_start      , ind_act_stop]      = regexp(str,'\$ACTU([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
[ind_vnav_start     , ind_vnav_stop]     = regexp(str,'\$VNMAR([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
if isempty(ind_vnav_start)
    [ind_vnav_start     , ind_vnav_stop] = regexp(str,'\$VNYPR([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
end
[ind_seg_start      , ind_seg_stop]      = regexp(str,'\$SEGM([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
[ind_spec_start     , ind_spec_stop]     = regexp(str,'\$SPEC([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
[ind_avgspec_start  , ind_avgspec_stop]  = regexp(str,'\$AVGS([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
[ind_dissrate_start , ind_dissrate_stop] = regexp(str,'\$RATE([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
[ind_apf0_start     , ind_apf0_stop]     = regexp(str,'\$APF0([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
[ind_apf1_start     , ind_apf1_stop]     = regexp(str,'\$APF1([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
[ind_apf2_start     , ind_apf2_stop]     = regexp(str,'\$APF2([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
[ind_fluor_start    , ind_fluor_stop]    = regexp(str,'\$ECOP([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
[ind_ttv_start      , ind_ttv_stop]      = regexp(str,'\$TTVP([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
if isempty(ind_ttv_start)
    [ind_ttv_start      , ind_ttv_stop]      = regexp(str,'\$TTV1([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
end
>>>>>>> Stashed changes

%$fluor0000000000045cf00000001c*6B0000000000045ceXXXXYYYYZZZZ*57
%% Define the header tag format
% In the versions of the SOM acquisition software since 23 May 2021 (and
% maybe before that?), these data types all have the same header tag
% format:
%     'EFE4','SB49','ALTI','ACTU','VNAV'

% get the offsets to parse the str
get_inds = @(x)(x.offset+1:x.offset+x.length);

tag.sync.strvalue='$';
tag.sync.offset=0;
tag.sync.length=1;
tag.sync.inds = get_inds(tag.sync);

tag.header.strvalue  = 'FFFF';
tag.header.length = strlength(tag.header.strvalue);
tag.header.offset = tag.sync.offset+tag.sync.length;
tag.header.inds = get_inds(tag.header);

tag.hextimestamp.strvalue  = "0000000000000000";
tag.hextimestamp.length = strlength(tag.hextimestamp.strvalue);
tag.hextimestamp.offset = tag.header.offset+tag.header.length;
tag.hextimestamp.inds = get_inds(tag.hextimestamp);

tag.hexlengthblock.strvalue  = "00000000";
tag.hexlengthblock.length = strlength(tag.hexlengthblock.strvalue);
tag.hexlengthblock.offset = tag.hextimestamp.offset+tag.hextimestamp.length;
tag.hexlengthblock.inds = get_inds(tag.hexlengthblock);

tag.headerchecksum.strvalue = "*FF";
tag.headerchecksum.length   = strlength(tag.headerchecksum.strvalue);
tag.headerchecksum.offset   = tag.hexlengthblock.offset+tag.hexlengthblock.length;
tag.headerchecksum.inds = get_inds(tag.headerchecksum);

tag.chksum.strvalue = "FFFFF" ;
tag.chksum.length   = strlength(tag.chksum.strvalue);

tag.data_offset = tag.headerchecksum.offset+tag.headerchecksum.length+1;

tag.laptoptime.strvalue = '0000000000';
tag.laptoptime.length = strlength(tag.laptoptime.strvalue);
tag.laptoptime.offset = -11;
tag.laptoptime.indx = get_inds(tag.laptoptime);

%% Start a list of processed and unprocessed data types
processed_data_types = {};
no_data_types = {};

%% Process setup SOM3 data
if isempty(ind_som_start)
    no_data_types = [no_data_types,'setup'];
    setup=[];
    setup=mod_som_read_setup_from_config('EPSILOMETER/config_files/bench_config');
    Meta_Data=epsiSetup_fill_meta_data(Meta_Data,setup);
else
    str_setup=str(ind_som_start+32:ind_som_stop-5);
    % str_setup=str(ind_som_start:ind_som_stop);
    setup=mod_som_read_setup_from_raw(str_setup);
    Meta_Data=epsiSetup_fill_meta_data(Meta_Data,setup);

    if ~isempty(epsi_probes)
        Meta_Data.AFE.t1.SN=epsi_probes.ch1.SN;
        Meta_Data.AFE.t1.cal=epsi_probes.ch1.cal;
        Meta_Data.AFE.t2.SN=epsi_probes.ch2.SN;
        Meta_Data.AFE.t2.cal=epsi_probes.ch2.cal;
        Meta_Data.AFE.s1.SN=epsi_probes.ch3.SN;
        Meta_Data.AFE.s2.SN=epsi_probes.ch4.SN;
        if epsi_probes.ch3.cal~=0
            Meta_Data.AFE.s1.cal=epsi_probes.ch3.cal;
        else
            AFE=get_shear_calibration(Meta_Data.AFE);
            Meta_Data.AFE=AFE;
            %ALB epsi case
            if isfield(Meta_Data.AFE.s1,'Sv')
                Meta_Data.AFE.s1.cal=Meta_Data.AFE.s1.Sv;
                Meta_Data.AFE.s2.cal=Meta_Data.AFE.s2.Sv;
            else
            %ALB fctd case
            % nothing to do so far 
            end
        end
    end
    
    if newSetup_flag
        Meta_Data.experiment_name    = experiment_name;
        Meta_Data.cruise_name        = cruise_name;
        Meta_Data.vehicle_name       = vehicle_name;
        Meta_Data.pressure_case_name = pressure_case_name;
        Meta_Data.fishflag_name      = fishflag_name;
    end

end
if isempty(ind_dcal_start)
end

%% GPS data We have to start with GPS because we need latitude 
if isempty(ind_gps_start)
    no_data_types = [no_data_types,'gps'];
    gps = [];
else
    processed_data_types = [processed_data_types,'gps'];
    %disp('processing gps data')
    
    % GPS-specific quantities
    % ---------------------------
    gpsmeta.data.n_blocks = numel(ind_gps_start);
    gpsmeta.data.recs_per_block = 1;
    gpsmeta.data.n_recs = gpsmeta.data.n_blocks*gpsmeta.data.recs_per_block;
    
    ind_time_start = ind_gps_start-10;
    ind_time_stop  = ind_gps_start-1;
    
    % San added these GPS steps. Here, he reads the file header to get system time.
    % You'll use this to correct the gps timestamp.
    FID = fopen(filename,'r');
    if FID<0
        error('MATLAB:FastCTD_ReadASCII:FileError', 'Could not open file %s',fname);
    end
    
    % NC 7/22/21 Problem! FastCTD_ASCII_parseheader doesn't work if there is
    % no header. It might be fixed by moving the "if there is gps data" line
    % above this step.
    FCTD = FastCTD_ASCII_parseheader(FID);
    
    %convert time to MATLAB time
    if FCTD.header.offset_time < 0
        FCTD.header.offset_time = correctNegativeTime(FCTD.header.offset_time)/86400+datenum(1970,1,1);
    else
        FCTD.header.offset_time = FCTD.header.offset_time/86400+datenum(1970,1,1);
    end
    if FCTD.header.system_time < 0
        FCTD.header.system_time = correctNegativeTime(FCTD.header.system_time)/86400/100+FCTD.header.offset_time;
    else
        FCTD.header.system_time = FCTD.header.system_time/86400/100+FCTD.header.offset_time;
    end
    
    fclose(FID);
    
    % Process GPS data
    % ---------------------------
    
    % Pre-allocate space for data
    gps.dnum   = nan(gpsmeta.data.n_recs,1);
    
    gps.latitude = nan(gpsmeta.data.n_recs,1);
    gps.longitude = nan(gpsmeta.data.n_recs,1);
    
    % Loop through data blocks and parse strings
    for iB=1:gpsmeta.data.n_blocks
        try
            % Grab the block of data starting with the header
            gps_block_str = str(ind_gps_start(iB):ind_gps_stop(iB)); %Moved here by Bethan June 26
            
            
            % Get the data after the header
            gps_block_data = str(ind_gps_start(iB):ind_gps_stop(iB));
            
            % Split the data into parts
            data_split = strsplit(gps_block_data,',');
            
            gps.dnum(iB) = str2double(str(ind_time_start(iB):ind_time_stop(iB)))/100/24/3600+FCTD.header.offset_time;
            
            gps.latitude(iB) = floor(str2double(data_split{3})/100) + mod(str2double(data_split{3}),100)/60; % then add minutes
            gps.latitude(iB) = gps.latitude(iB).*(2*strcmpi(data_split{4},'N')-1); % check for north or south
            
            gps.longitude(iB) = floor(str2double(data_split{5})/100) + mod(str2double(data_split{5}),100)/60; % then add minutes
            gps.longitude(iB) = gps.longitude(iB).*(2*strcmpi(data_split{6},'E')-1); % check for East or West (if west multiply by -1)
        catch err
            iB
        end
    end
    
end %end loop if there is gps data


%% Process EFE data
if isempty(ind_efe_start)
    no_data_types = [no_data_types,'epsi'];
    epsi=[];
else
    processed_data_types = [processed_data_types,'epsi'];
    %disp('processing epsi data')
    % EFE-specific quantities
    efe.data.n_blocks           = numel(ind_efe_start); %Number of data blocks beginning with $EFE
    
    efe.data.n_channels         = Meta_Data.PROCESS.nb_channels;
    if efe.data.n_channels==3
        efe.data.sample_freq        = 160;
    else
        efe.data.sample_freq        = 320;
    end
    efe.data.sample_period      = 1/efe.data.sample_freq;
    efe.data.bytes_per_channel  = 3;
    efe.data.timestamp_length   = 8;
    efe.data.n_elements         = efe.data.timestamp_length + ...
        efe.data.n_channels*efe.data.bytes_per_channel; % 8 bytes timestamps + 3 bytes ADC
    efe.data.recs_per_block     = 80;
    efe.data.n_recs             = efe.data.n_blocks*efe.data.recs_per_block;
    
    efe.bit_counts              = 24;
    efe.gain                    = 1;
    efe.acc_offset              = 1.8/2;
    % NC on 16 July 2021 during BLT cruise - changed acc_factor from 0.5 to
    % 0.4. The accelerometer gain from the manual is 400 mV/g
    %efe.acc_factor              = 0.5;
    efe.acc_factor              = 0.4;
    efe.channels                = Meta_Data.PROCESS.channels;
    
    % These fields appear to be unused:
    %efe.data.nblocks            = 160;
    %efe.data.length = efe.data.nblocks*efe.data.nchannels*efe.data.bytes_per_channel;
    %efe.data.n_recs = numel(ind_efe_start);
    
    
    % Pre-allocate space for data
    epsi_timestamp   = nan(efe.data.n_recs,1);
    %epsi.time_s and epsi.dnum will be created from epsi_timestamp once
    %all its records are filled
    
    % Loop through data blocks and parse strings
    for iB = 1:efe.data.n_blocks
        
        % Grab the block of data starting with the header
        efe_block_str = str(ind_efe_start(iB):ind_efe_stop(iB));
        
        % Get the header timestamp and length of data block (?)
        efe.data.hextimestamp.value   = hex2dec(efe_block_str(tag.hextimestamp.inds));
        efe.data.hexlengthblock.value = hex2dec(efe_block_str(tag.hexlengthblock.inds));
        % It should be the case that
        % efe.hexlengthblock.value==efe.data.element_length*efe.data.recs_per_block
        
        % Get the data after the header.
        efe_block_data = efe_block_str(tag.data_offset:end-tag.chksum.length);
        
        % Does length(efe_block_data)=efe.hexlengthblock.value?
        if (length(efe_block_data)~=efe.data.n_elements*efe.data.recs_per_block)
            fprintf("EFE block %i has incorrect length\r\n",iB)
        else
            
            % Reshape the data
            efe.data.raw_bytes{iB} = reshape(uint32(efe_block_data),...
                efe.data.n_elements,...
                efe.data.recs_per_block).';
            
            % Get the checksum value
            i1 = tag.data_offset+efe.data.hexlengthblock.value;
            i2 = i1+tag.chksum.length-1;
            % The full checksum looks like *8F__ so we get str indices from
            % i1+1 to i2-2 to isolate just the 8F part
            efe.checksum.data(iB) = hex2dec(efe_block_str(i1+1:i2-2));
            
        end %end if efe data block is the correct size
    end %end loop through efe blocks
    
    if ~isfield(efe.data,'raw_bytes')
        no_data_types = [no_data_types,'epsi'];
        epsi=[];
    elseif isfield(efe.data,'raw_bytes')
        % Clean and concatenate efe.data.raw_bytes
        % Check for empty raw_bytes cells and concatenate only the full ones
        idnull=cellfun(@isempty,efe.data.raw_bytes);
        efe.data.raw_bytes=cell2mat(efe.data.raw_bytes(~idnull).');
        
        % The first 8 bytes are the timestamp. Timestamp in milliseconds since
        % 1970 is too big for uint32, so use uint64
        ind_time = 1:efe.data.timestamp_length;
        mult = 256.^[0:7].';
        epsi_timestamp  = double(uint64(efe.data.raw_bytes(:,ind_time)))*mult;
        
        % Get ADC data. After the timestamp, data are in 3-byte units
        byte1 = efe.data.raw_bytes(:,ind_time(end)+1:3:end);
        byte2 = efe.data.raw_bytes(:,ind_time(end)+2:3:end);
        byte3 = efe.data.raw_bytes(:,ind_time(end)+3:3:end);
        efe.data.raw = byte1*256^2 + byte2*256^1 + byte3*256^0;
        
        
        % Save data as 'counts'
        for cha=1:efe.data.n_channels
            wh_channel=efe.channels{cha};
            epsi.([wh_channel '_count']) = efe.data.raw(:,cha);
        end
        
        % Convert counts to volts (for s1,s2,t1,t2) and g-units (for a1,a2,a3)
        Unipolar=@(FR,data) (FR/efe.gain*double(data)/2.^(efe.bit_counts));
        Bipolar=@(FR,data) (FR/efe.gain*(double(data)/2.^(efe.bit_counts-1)-1));
        %     % NC 17 July 2021 - the Bipolar function was mapping data to -2.5 to 0
        %     % V. We want to map it to -1.25 to 1.25 V. JUST KIDDING! We convinced
        %     ourselves  it should indeed by to -2.5 to 2.5 V
        %     Bipolar = @(FR,data) (FR/efe.gain/2*(double(data)/2.^(efe.bit_counts-1)-1));
        
        for cha=1:efe.data.n_channels
            wh_channel=efe.channels{cha};
            FR=Meta_Data.AFE.(wh_channel).full_range;
            switch Meta_Data.AFE.(wh_channel).ADCconf
                case {'Bipolar','bipolar'}
                    epsi.([wh_channel '_volt'])=Bipolar(FR,epsi.([wh_channel '_count']));
                case {'Unipolar','unipolar'}
                    epsi.([wh_channel '_volt'])=Unipolar(FR,epsi.([wh_channel '_count']));
                    
            end
            
            switch Meta_Data.AFE.(wh_channel).full_range
                case 1.8 %case for acceleration channels
                    epsi.([wh_channel '_g']) = (epsi.([wh_channel '_volt'])-efe.acc_offset)/efe.acc_factor;
                    epsi = rmfield(epsi,[wh_channel '_volt']);
            end
        end
        
        % Reshape to column arrays - must have been left from a previous
        % version. They are already column arrays
        %     epsiFields = fieldnames(epsi);
        %     for iF  = 1:numel(epsiFields)
        %         epsi.(epsiFields{iF}) = reshape(epsi.(epsiFields{iF}),[],1);
        %     end
        
        % If timestamp has values like 1.6e12, it is in milliseconds since Jan
        % 1, 1970. Otherwise it's in milliseconds since the start of the record
        if nanmedian(epsi_timestamp)>1e9
            % time_s - seconds since 1970
            % dnum - matlab datenum
            [epsi.time_s,epsi.dnum] = convert_timestamp(epsi_timestamp);
        else
            % time_s - seconds since power on
            epsi.time_s = epsi_timestamp./1000;
            epsi.dnum = Meta_Data.start_dnum + days(seconds(epsi.time_s));
        end
        
        % Sort epsi fields
        if contains(Meta_Data.fishflag_name,'FCTD')
            epsi = orderfields(epsi,{'dnum','time_s','t1_count','t2_count','f1_count',...
                'c1_count','a1_count','a2_count','a3_count','t1_volt','t2_volt','f1_volt',...
                'c1_volt','a1_g','a2_g','a3_g'});
        else
        try
            epsi = orderfields(epsi,{'dnum','time_s','t1_count','t2_count','s1_count',...
                's2_count','a1_count','a2_count','a3_count','t1_volt','t2_volt','s1_volt',...
                's2_volt','a1_g','a2_g','a3_g'});
        catch
            epsi = orderfields(epsi,{'dnum','time_s','t1_count','s1_count',...
                'a3_count','a2_count','t1_volt','s1_volt',...
                'a3_g','a2_g'});
        end
        end%end of contains FCTD
    end %end if isfield(efe.data,'raw_bytes')
end %end if there is epsi data


%% CTD
if isempty(ind_sbe_start)
    no_data_types = [no_data_types,'ctd'];
    ctd=[];
else
    processed_data_types = [processed_data_types,'ctd'];
    bad_SB41sample_flag=0;
    %disp('processing ctd data')
    
    % SBE-specific quantities
    % ---------------------------------
    
    switch Meta_Data.CTD.name
        case{"SBE49","SBE","S49","SB49"}
            sbe.data.format      = 'eng';
            sbe.data.length      = 24; %It's 24 for BLT data
            sbe.data.sample_freq = 16;
            sbe.cal              = SBEcal; % SBEcal is readfrom the Header of the modraw
        case{"SBE41","S41","SB41"}
            sbe.data.format      = 'PTS';
            sbe.data.length      = 28;
            sbe.data.sample_freq = 1;
            % sbe.cal = ? Is calibration calculated differently for these
            % instruments?
    end
    
    sbe.data.sample_period    = 1/sbe.data.sample_freq;
    sbe.data.n_blocks         = numel(ind_sbe_start);
    sbe.data.recs_per_block   = Meta_Data.CTD.sample_per_record;
    sbe.data.n_recs           = sbe.data.n_blocks*sbe.data.recs_per_block;
    sbe.data.timestamp_length = 16;
    
    % These fields appear to be unused:
    %sbe.block_time            = nan(sbe.data.n_blocks,1);
    %sbe.checksum              = uint8(zeros(sbe.data.n_recs,1));
    
    
    
    
    % Process SBE data
    % --------------------------------
    
    % Pre-allocate space for data
    ctd_timestamp   = nan(sbe.data.n_recs,1);
    %ctd.ctdtime and ctd.ctddnum will be created from ctd_timestamp once
    %all its records are filled
    ctd.P_raw       = nan(sbe.data.n_recs,1);
    ctd.T_raw       = nan(sbe.data.n_recs,1);
    ctd.S_raw       = nan(sbe.data.n_recs,1);
    ctd.C_raw       = nan(sbe.data.n_recs,1);
    ctd.PT_raw      = nan(sbe.data.n_recs,1);
    ctd.P           = nan(sbe.data.n_recs,1);
    ctd.T           = nan(sbe.data.n_recs,1);
    ctd.S           = nan(sbe.data.n_recs,1);
    ctd.C           = nan(sbe.data.n_recs,1);
    
    % Initialize datarecord counter
    n_rec = 0;
    
    % Loop through data blocks and parse strings
    for iB = 1:sbe.data.n_blocks
        
        % Grab the block of data starting with the header
        sbe_block_str = str(ind_sbe_start(iB):ind_sbe_stop(iB));
        
        % Get the header timestamp and length of data block
        sbe.hextimestamp.value   = hex2dec(sbe_block_str(tag.hextimestamp.inds));
        sbe.hexlengthblock.value = hex2dec(sbe_block_str(tag.hexlengthblock.inds));
        % It should be the case that
        % sbe.hexlengthblock.value==(24+16)*sbe.data.recs_per_block
        % Where do 16 and 24 come from?
        
        % Get the data after the header.
        sbe_block_data = sbe_block_str(tag.data_offset:end-tag.chksum.length);
        
        % Where do 16 and 24 come from? % Does length(sbe_block_data)=sbe.hexlengthblock.value?
        if (length(sbe_block_data)~=(sbe.data.length+16)*sbe.data.recs_per_block)
            fprintf("SBE block %i has incorrect length\r\n",iB)
        else
            %
            % Where do 16 and 24 come from?
            sbe_block_data=reshape(sbe_block_data,16+sbe.data.length,sbe.data.recs_per_block).';
            
            for iR=1:sbe.data.recs_per_block
                
                % Count up the record number
                n_rec=n_rec+1;
                
                % The ctd data is everything but the last two elements of the
                % block
                element_ctd=sbe_block_data(iR,1:end-2);
                
                % The hexadecimal timestamp is the first 16 characters of the
                % data
                ctd_timestamp(n_rec) = hex2dec(element_ctd(1:16));
                
                % Everything after that is the data
                rec_ctd=element_ctd(sbe.data.timestamp_length+1:end);
                
                % Split the data into parts and convert according to
                % sbe.data.format
                switch sbe.data.format
                    case 'PTS'
                        data_split        = strsplit(rec_ctd,',');
                        if length(data_split)==3
                            ctd.P(n_rec)      = str2double(data_split{1});
                            ctd.T(n_rec)      = str2double(data_split{2});
                            ctd.S(n_rec)      = str2double(data_split{3});
                            ctd.C(n_rec)      = NaN;
                        else
                            bad_SB41sample_flag=1;
                        end
                    case 'eng'
                        ctd.T_raw(n_rec)  = hex2dec(rec_ctd(:,1:6));
                        ctd.C_raw(n_rec)  = hex2dec(rec_ctd(:,(1:6)+6));
                        ctd.P_raw(n_rec)  = hex2dec(rec_ctd(:,(1:6)+12));
                        ctd.PT_raw(n_rec) = hex2dec(rec_ctd(:,(1:4)+18));
                end %end switch sbe.data.format
            end %end loop through recs_per_block
            
            % If timestamp has values like 1.6e12, it is in milliseconds since Jan
            % 1, 1970. Otherwise it's in milliseconds since the start of the record
            if nanmedian(ctd_timestamp)>1e9
                % time_s - seconds since 1970
                % dnum - matlab datenum
                [ctd.time_s,ctd.dnum] = convert_timestamp(ctd_timestamp);
            else
                % time_s - seconds since power on
                ctd.time_s = ctd_timestamp./1000;
                ctd.dnum = Meta_Data.start_dnum + days(seconds(ctd.time_s));
            end
            
            % If the data were in engineering units, convert to physical units
            switch sbe.data.format
                case 'eng'
                    ctd = sbe49_ascii_get_temperature(ctd,sbe);
                    ctd = sbe49_ascii_get_pressure(ctd,sbe);
                    ctd = sbe49_ascii_get_conductivity(ctd,sbe);
                    
                    ctd.S    = sw_salt(ctd.C*10./c3515,ctd.T,ctd.P);
            end
            
            ctd.th   = sw_ptmp(ctd.S,ctd.T,ctd.P,0);
            ctd.sgth  = sw_pden(ctd.S,ctd.T,ctd.P,0);
            ctd.dPdt = [0; diff(ctd.P)./diff(ctd.time_s)];
            
            % NC 17 July 2021 - added ctd.z and ctd.dzdt.
            % get_scan_spectra.m will use dzdt to define fall speed w.
            if ~isfield(Meta_Data.PROCESS,'latitude')
                error('Need latitude to get depth from pressure data. Add to MetaProcess text file.')
            else
                % ALB add the gps data if it exists. 
                if isfield(gps,'latitude')
                    try
                        if length(gps.latitude)>1
                            ctd_lat  = interp1(gps.dnum,gps.latitude,ctd.dnum);
                            ctd.z    = sw_dpth(ctd.P,ctd_lat);
                        else
                            ctd_lat  = gps.latitude;
                            ctd.z    = sw_dpth(ctd.P,ctd_lat);
                        end
                    catch
                        warning("Can not interpolate latitude ctd.num. Line 665 mod_som_read_epsi_v4")
                    end
                else
                    ctd.z    = sw_dpth(ctd.P,Meta_Data.PROCESS.latitude);
                end
                try
                ctd.dzdt = [0; diff(ctd.z)./diff(ctd.time_s)];
                catch
                    disp('issue with dzdt mod_som_read_epsi_files_v4 l.689')
                end
            end
            
            if bad_SB41sample_flag==0
                % Make sure ctd.S is real. Every once in a while, S comes out imaginary, I think only when SBE is on deck.
                ctd.S = real(ctd.S);
                
                % Sort ctd fields
                ctd = orderfields(ctd,{'dnum','time_s','P_raw','T_raw','S_raw',...
                    'C_raw','PT_raw','P','z','T','S','C','th','sgth','dPdt','dzdt'});
            end
        end %end if sbe data block is the correct size
%         % Tracking blocks
%         if mod(n_rec,1000)==0
%             sprintf("block %i \r\n",n_rec)
%         end
    end %end loop through sbe blocks
    
    
end %end loop if there is ctd data

%% MOD  Altimeter
if isempty(ind_alt_start)
    no_data_types = [no_data_types,'alt'];
    alt=[];
else
    processed_data_types = [processed_data_types,'alt'];
    %disp('processing alt data')
    
    % ALTI-specific quantities
    % --------------------------------
    alti.data.n_blocks       = numel(ind_alt_start);
    alti.data.recs_per_block = 1;
    alti.data.n_recs         = alti.data.n_blocks*alti.data.recs_per_block;
    
    alti.ten_microsec2sec    = 1e-5;
    alti.sound_speed         = 1500;
    
    % Process ALTI data
    % --------------------------------
    
    % Pre-allocate space for data
    alt_timestamp  = nan(alti.data.n_recs,1);
    %alt.time_s and alt.dnum will be created from alt_timestamp once
    %all its records are filled
    alt.dst        = nan(alti.data.n_recs,1);
    
    % Loop through data blocks and parse strings
    for iB=1:alti.data.n_blocks
        
        % Grab the block of data starting with the header
        alt_block_str = str(ind_alt_start(iB):ind_alt_stop(iB));
        
        % For the altimeter, all the data is actually in the header
        try
            alti.hextimestamp.value   = hex2dec(alt_block_str(tag.hextimestamp.inds));
        catch
            alti.hextimestamp.value = nan;
        end
        try
            alt_timestamp(iB) = alti.hextimestamp.value;
        catch
            alt_timestamp(iB) = nan;
        end
        
        % The altimeter block does not have a hexlengthblock (hexadecimal
        % length of data block). It has the altimeter reading in that
        % position
        alti.hexlengthblock.value = alt_block_str(tag.hexlengthblock.inds);
        dst_microsec = str2double(alti.hexlengthblock.value);
        alt.dst(iB) = dst_microsec*alti.ten_microsec2sec*alti.sound_speed;
        
        % Calculate actual height above bottom (hab) from altimeter
        % distance reading. The altimeter is angled at 10 degrees and is positioned 5 ft above the
        % crashguard. The probes sit 2.02 inches behind the crash guard.
        % (See Epsi Processing Manual - Altimeter correction for diagram).
        convert_dissrate = @(x) ((x-apf.data.dissrate_count0)/apf.data.dissrate_per_bit);
        feet2meters = @(x) (x*0.3048);
        inches2meters = @(x) (x*0.0254);
        angle_deg = Meta_Data.PROCESS.alt_angle_deg;
        alt_to_crashguard = Meta_Data.PROCESS.alt_dist_from_crashguard_ft;
        probe_to_crashguard = Meta_Data.PROCESS.alt_probe_dist_from_crashguard_in;

        A = alt.dst(iB);
        theta = deg2rad(angle_deg);
        altimeter_height_above_probes = feet2meters(alt_to_crashguard) - ...
            inches2meters(probe_to_crashguard);

        % alt.hab is the height of the probes above the bottom
        alt.hab(iB,1) = A*cos(theta) - altimeter_height_above_probes;

    end %end loop through alt blocks
    
    % If timestamp has values like 1.6e12, it is in milliseconds since Jan
    % 1, 1970. Otherwise it's in milliseconds since the start of the record
    if nanmedian(alt_timestamp)>1e9
        % time_s - seconds since 1970
        % dnum - matlab datenum
        [alt.time_s,alt.dnum] = convert_timestamp(alt_timestamp);
    else
        % time_s - seconds since power on
        alt.time_s = alt_timestamp./1000;
        alt.dnum = Meta_Data.start_dnum + days(seconds(alt.time_s));
    end
    
    % Order alt fields
    alt = orderfields(alt,{'dnum','time_s','dst','hab'});
    
end %end loop if there is alt data

%% ISA500
if isempty(ind_isap_start)
    no_data_types = [no_data_types,'isap'];
    isap=[];
else
    processed_data_types = [processed_data_types,'isap'];
    %disp('processing alt data')
    
    % ALTI-specific quantities
    % --------------------------------
    isa.data.n_blocks       = numel(ind_isap_start);
    isa.data.recs_per_block = 1;
    isa.data.n_recs         = isa.data.n_blocks*isa.data.recs_per_block;
    
    isa.ten_microsec2sec    = 1e-5;
    isa.sound_speed         = 1500;
    
    % Process ISAP data
    % --------------------------------
    
    % Pre-allocate space for data
    isap_timestamp  = nan(isa.data.n_recs,1);
    %alt.time_s and alt.dnum will be created from alt_timestamp once
    %all its records are filled
    isap.dst        = nan(isa.data.n_recs,1);
    
    % Loop through data blocks and parse strings
    for iB=1:isa.data.n_blocks
        
        % Grab the block of data starting with the header
        isap_block_str = str(ind_isap_start(iB):ind_isap_stop(iB));
        
        % For the altimeter, all the data is actually in the header
        try
            isa.hextimestamp.value   = hex2dec(isap_block_str(tag.hextimestamp.inds));
        catch
            isa.hextimestamp.value = nan;
        end
        try
            isap_timestamp(iB) = isa.hextimestamp.value;
        catch
            isap_timestamp(iB) = nan;
        end
        
        % The altimeter block does not have a hexlengthblock (hexadecimal
        % length of data block). It has the altimeter reading in that
        % position
        isa.hexlengthblock.value = isap_block_str(tag.hexlengthblock.inds);
        
        isap.dst(iB) = str2double(isap_block_str(tag.data_offset+ ...
                                 tag.hextimestamp.length+ ...
                                 tag.header.offset+ ...
                                 tag.header.length+2:end-5));
        % ALB ISAP500 send meters directly
        % isap.dst(iB) = dst_microsec*isap.ten_microsec2sec*isap.sound_speed;
        
        % Calculate actual height above bottom (hab) from altimeter
        % distance reading. The altimeter is angled at 10 degrees and is positioned 5 ft above the
        % crashguard. The probes sit 2.02 inches behind the crash guard.
        % (See Epsi Processing Manual - Altimeter correction for diagram).
        % convert_dissrate = @(x) ((x-apf.data.dissrate_count0)/apf.data.dissrate_per_bit);
        feet2meters = @(x) (x*0.3048);
        inches2meters = @(x) (x*0.0254);
        angle_deg = Meta_Data.PROCESS.alt_angle_deg;
        isap_to_crashguard = Meta_Data.PROCESS.alt_dist_from_crashguard_ft;
        probe_to_crashguard = Meta_Data.PROCESS.alt_probe_dist_from_crashguard_in;

        A = isap.dst(iB);
        theta = deg2rad(angle_deg);
        altimeter_height_above_probes = feet2meters(isap_to_crashguard) - ...
            inches2meters(probe_to_crashguard);

        % alt.hab is the height of the probes above the bottom
        isap.hab(iB,1) = A*cos(theta) - altimeter_height_above_probes;

    end %end loop through alt blocks
    
    % If timestamp has values like 1.6e12, it is in milliseconds since Jan
    % 1, 1970. Otherwise it's in milliseconds since the start of the record
    if nanmedian(isap_timestamp)>1e9
        % time_s - seconds since 1970
        % dnum - matlab datenum
        [isap.time_s,isap.dnum] = convert_timestamp(isap_timestamp);
    else
        % time_s - seconds since power on
        isap.time_s = isap_timestamp./1000;
        isap.dnum = Meta_Data.start_dnum + days(seconds(isap.time_s));
    end
    
    % Order alt fields
    isap = orderfields(isap,{'dnum','time_s','dst','hab'});
    
end %end loop if there is alt data


%% Actuator
if isempty(ind_act_start)
    no_data_types = [no_data_types,'act'];
    act=[];
else
    processed_data_types = [processed_data_types,'act'];
    %disp('processing act data')
    
    % ACTU-specific quantities
    % ---------------------------
    
    
    
    % Process ACTU data
    % ---------------------------
    
    % Pre-allocate space for data
    %   act_timestamp   = nan(actu.data.n_recs,1);
    %act.acttime and act.actdnum will be created from act_timestamp once
    %all its records are filled
    
    act = [];
end

%% Vector navigation

if isempty(ind_vnav_start)
    no_data_types = [no_data_types,'vnav'];
    vnav=[];
else
    processed_data_types = [processed_data_types,'vnav'];
    %disp('processing vnav data')
    
    % VNAV-specific quantities
    % ---------------------------
    vecnav.data.n_blocks = numel(ind_vnav_start);
    vecnav.data.recs_per_block = 1;
    vecnav.data.n_recs = vecnav.data.n_blocks*vecnav.data.recs_per_block;
    
    ind_time_start = ind_vnav_start-16;
    ind_time_stop  = ind_vnav_start-1;
    
    % Process VNAV data
    % ---------------------------
    % Pre-allocate space for data
    vnav_timestamp   = nan(vecnav.data.n_recs,1);
    %vnav.time_s and vnav.dnum will be created from vnav_timestamp once
    %all its records are filled
    vnav.compass = nan(vecnav.data.n_recs,3);
    vnav.acceleration = nan(vecnav.data.n_recs,3);
    vnav.gyro = nan(vecnav.data.n_recs,3);
    vnav.yaw = nan(vecnav.data.n_recs,1);
    vnav.pitch = nan(vecnav.data.n_recs,1);
    vnav.roll = nan(vecnav.data.n_recs,1);
    
    %     % Grab the block of data starting with the header
    %     vnav_block_str = str(ind_vnav_start(iB):ind_vnav_stop(iB));
    %     %Commented out and moved below by Bethan June 26
    Lsample=length(str(ind_vnav_start(iB):ind_vnav_stop(iB)));
    
    % Loop through data blocks and parse strings
    for iB=1:vecnav.data.n_blocks
        
        % Grab the block of data starting with the header
        vnav_block_str = str(ind_vnav_start(iB):ind_vnav_stop(iB)); %Moved here by Bethan June 26
        % For the vecnav, we use the hexadecimal timestamp that comes before
        % $VNMAR
        % vnav_timestamp(iB) = hex2dec(str(ind_time_start(iB):ind_time_stop(iB)));
        vnav_timestamp(iB) = hex2dec(str(ind_vnav_start(iB)-16:ind_vnav_start(iB)-1));
        
        % Get the data after the header
        vnav_block_data = str(ind_vnav_start(iB):ind_vnav_stop(iB)-tag.chksum.length);
        
        % Split the data into parts
        data_split = strsplit(vnav_block_data,',');
        
        % Compass, acceleration, and gyro each have x, y, and z components
        switch Lsample
            case 38 %yaw pitch roll
                vnav.yaw(iB)   = str2double(data_split{2});
                vnav.pitch(iB) = str2double(data_split{3});
                vnav.roll(iB)  = str2double(data_split{4});
            otherwise
                for iC=1:3
                    vnav.compass(iB,iC)=str2double(data_split{1,iC+1}); %Compass (x,y,z) units = gauss
                    vnav.acceleration(iB,iC)=str2double(data_split{1,iC+4}); %Acceleration (x,y,x) units = m/s^2
                    vnav.gyro(iB,iC)=str2double(data_split{1,iC+7}); %Gyro (x,y,z) units = rad/s
                end
        end
        
    end
    
    % If timestamp has values like 1.6e12, it is in milliseconds since Jan
    % 1, 1970. Otherwise it's in milliseconds since the start of the record
    if nanmedian(vnav_timestamp)>1e9
        % time_s - seconds since 1970
        % dnum - matlab datenum
        [vnav.time_s,vnav.dnum] = convert_timestamp(vnav_timestamp);
    else
        % time_s - seconds since power on
        vnav.time_s = vnav_timestamp./1000;
        vnav.dnum = Meta_Data.start_dnum + days(seconds(vnav.time_s));
    end
    
    % Order vnav fields
    vnav = orderfields(vnav,{'dnum','time_s','compass','acceleration','gyro','yaw','pitch','roll'});
    
end %end loop if there is vnav data


%% Process SEGM data
if isempty(ind_seg_start)
    no_data_types = [no_data_types,'seg'];
    seg=[];
else
    processed_data_types = [processed_data_types,'seg'];
    %disp('processing segment data')
    
    % Pre-allocate space for data
    seg.data.n_blocks           = numel(ind_seg_start); %Number of data blocks beginning with $EFE
    
    seg.data.n_channels         = 3;
    seg.data.sample_freq        = 320;
    seg.data.sample_period      = 1/seg.data.sample_freq;
    seg.data.sample_per_segment = 2048; %NFFT
    seg.data.segment_per_block  = 1;
    seg.data.timestamp_length   = 8;
    seg.data.n_elements         = seg.data.timestamp_length+ ...
        seg.data.n_channels*...
        seg.data.sample_per_segment*4;
    seg.channels ={'t1_volt','s1_volt','a3_g'};
    
    %     seg.data.n_recs             = efe.data.n_blocks*efe.data.recs_per_block;
    % The first 8 bytes are the timestamp. Timestamp in milliseconds since
    % 1970 is too big for uint32, so use uint64
    ind_time = 1:seg.data.timestamp_length;
    mult = 256.^[0:7].';
    
    
    segment_timestamp           = nan(seg.data.n_blocks,1);
    %epsi.time_s and epsi.dnum will be created from epsi_timestamp once
    %all its records are filled
    
    % Loop through data blocks and parse strings
    for iB = 1:seg.data.n_blocks
        
        % Grab the block of data starting with the header
        seg_block_str = str(ind_seg_start(iB):ind_seg_stop(iB));
        
        % Get the header timestamp and length of data block (?)
        seg.data.hextimestamp.value   = hex2dec(seg_block_str(tag.hextimestamp.inds));
        seg.data.hexlengthblock.value = hex2dec(seg_block_str(tag.hexlengthblock.inds));
        % It should be the case that
        % efe.hexlengthblock.value==efe.data.element_length*efe.data.recs_per_block
        
        % Get the data after the header.
        seg_block_data = seg_block_str(tag.data_offset:end-tag.chksum.length);
        
        % Does length(efe_block_data)=efe.hexlengthblock.value?
        if (length(seg_block_data)~=seg.data.n_elements)
            fprintf("SEG block %i has incorrect length\r\n",iB)
        else
            
            segment_timestamp(iB)=double(uint64(seg_block_data(ind_time)))*mult;
            
            % Reshape the data
            seg.data.raw_bytes{iB} = reshape(uint8(seg_block_data(ind_time(end)+1:end)),...
                4,...
                seg.data.sample_per_segment*...
                seg.data.n_channels);
            
            for n=1:seg.data.sample_per_segment*seg.data.n_channels
                seg.data.raw_float{iB}(n) = double(typecast(seg.data.raw_bytes{iB}(:,n).' , 'single'));
                %                 seg.data.raw_float{iB}(n) = double(uint32(seg.data.raw_bytes{iB}(:,n))).'*mult(1:4);
            end
            % Save data as 'counts'
            for cha=1:seg.data.n_channels
                wh_channel=seg.channels{cha};
                seg.(wh_channel){iB} = seg.data.raw_float{iB}(1+(cha-1)*seg.data.sample_per_segment:cha*seg.data.sample_per_segment);
            end
            
            % Get the checksum value
            i1 = tag.data_offset+seg.data.hexlengthblock.value;
            i2 = i1+tag.chksum.length-1;
            % The full checksum looks like *8F__ so we get str indices from
            % i1+1 to i2-2 to isolate just the 8F part
            seg.checksum.data(iB) = hex2dec(seg_block_str(i1+1:i2-2));
        end %end if efe data block is the correct size
    end %end loop through efe blocks
    
    % Clean and concatenate efe.data.raw_bytes
    % Check for empty raw_bytes cells and concatenate only the full ones
    for cha=1:seg.data.n_channels
        wh_channel=seg.channels{cha};
        idnull=cellfun(@isempty,seg.(wh_channel));
        %         seg.(wh_channel)=cell2mat(seg.(wh_channel)(~idnull));
        seg.(wh_channel)=seg.(wh_channel)(~idnull);
    end
    
    if nanmedian(segment_timestamp(iB))>1e9
        % time_s - seconds since 1970
        % dnum - matlab datenum
        [seg.time_s,seg.dnum] = convert_timestamp(segment_timestamp);
    else
        % time_s - seconds since power on
        seg.time_s = arrayfun(@(x,y) linspace(x,x+1000*(seg.data.sample_per_segment+1)/seg.data.sample_freq ...
            ,seg.data.sample_per_segment)./1000,segment_timestamp,'un',0);
        
        %         seg.time_s=cell2mat(seg.time_s.');
        seg.dnum = cellfun(@(x) Meta_Data.start_dnum + days(seconds(x)),seg.time_s,'un',0);
    end
    
    % Sort epsi fields
    seg = orderfields(seg,{'dnum','time_s','t1_volt','s1_volt','a3_g','data','channels','checksum'});
    
end %end if there is epsi data

%% Process SPEC data
if isempty(ind_spec_start)
    no_data_types = [no_data_types,'spec'];
    spec=[];
else
    try
    processed_data_types = [processed_data_types,'spec'];
    %disp('processing spec data')
    
    % Pre-allocate space for data
    spec.data.n_blocks           = numel(ind_spec_start); %Number of data blocks beginning with $EFE
    
    spec.data.n_channels          = 3;
    spec.data.sample_freq       = 320;
    spec.data.sample_per_spec   = 2048/2; %NFFT/2
    spec.data.spec_per_block  = 1;
    spec.data.timestamp_length   = 8;
    spec.data.n_elements         = spec.data.timestamp_length+ ...
        spec.data.n_channels*...
        spec.data.sample_per_spec*4;
    spec.channels ={'t1_volt','s1_volt','a3_g'};
    
    %     seg.data.n_recs             = efe.data.n_blocks*efe.data.recs_per_block;
    % The first 8 bytes are the timestamp. Timestamp in milliseconds since
    % 1970 is too big for uint32, so use uint64
    ind_time = 1:spec.data.timestamp_length;
    mult = 256.^[0:7].';
    
    
    spec_timestamp           = nan(spec.data.n_blocks,1);
    %epsi.time_s and epsi.dnum will be created from epsi_timestamp once
    %all its records are filled
    
    % Loop through data blocks and parse strings
    for iB = 1:spec.data.n_blocks
        
        % Grab the block of data starting with the header
        spec_block_str = str(ind_spec_start(iB):ind_spec_stop(iB));
        
        % Get the header timestamp and length of data block (?)
        spec.data.hextimestamp.value   = hex2dec(spec_block_str(tag.hextimestamp.inds));
        spec.data.hexlengthblock.value = hex2dec(spec_block_str(tag.hexlengthblock.inds));
        % It should be the case that
        % efe.hexlengthblock.value==efe.data.element_length*efe.data.recs_per_block
        
        % Get the data after the header.
        spec_block_data = spec_block_str(tag.data_offset:end-tag.chksum.length);
        
        % Does length(efe_block_data)=efe.hexlengthblock.value?
        if (length(spec_block_data)~=spec.data.n_elements)
            fprintf("SPEC block %i has incorrect length\r\n",iB)
        else
            
            spec_timestamp(iB)=double(uint64(spec_block_data(ind_time)))*mult;
            
            % Reshape the data
            spec.data.raw_bytes{iB} = reshape(uint8(spec_block_data(ind_time(end)+1:end)),...
                4,...
                spec.data.sample_per_spec*...
                spec.data.n_channels);
            
            for n=1:spec.data.sample_per_spec*spec.data.n_channels
                spec.data.raw_float{iB}(n) = double(typecast(spec.data.raw_bytes{iB}(:,n).' , 'single'));
                %                 seg.data.raw_float{iB}(n) = double(uint32(seg.data.raw_bytes{iB}(:,n))).'*mult(1:4);
            end
            % Save data as 'counts'
            for cha=1:spec.data.n_channels
                wh_channel=spec.channels{cha};
                spec.(wh_channel){iB} = spec.data.raw_float{iB}(1+(cha-1)*spec.data.sample_per_spec:cha*spec.data.sample_per_spec);
            end
            
            % Get the checksum value
            i1 = tag.data_offset+spec.data.hexlengthblock.value;
            i2 = i1+tag.chksum.length-1;
            % The full checksum looks like *8F__ so we get str indices from
            % i1+1 to i2-2 to isolate just the 8F part
            spec.checksum.data(iB) = hex2dec(spec_block_str(i1+1:i2-2));
        end %end if efe data block is the correct size
    end %end loop through efe blocks
    
    % Clean and concatenate efe.data.raw_bytes
    % Check for empty raw_bytes cells and concatenate only the full ones
    for cha=1:spec.data.n_channels
        wh_channel=spec.channels{cha};
        idnull=cellfun(@isempty,spec.(wh_channel));
        spec.(wh_channel)=spec.(wh_channel)(~idnull);
    end
    [~,fe]=pwelch(1.*(1:Meta_Data.PROCESS.nfft),Meta_Data.PROCESS.nfft,[], ...
        Meta_Data.PROCESS.nfft,Meta_Data.AFE.FS);
    spec.freq=fe(2:end);

    
    if nanmedian(spec_timestamp)>1e9
        % time_s - seconds since 1970
        % dnum - matlab datenum
        [spec.time_s,spec.dnum] = convert_timestamp(spec_timestamp);
    else
        % time_s - seconds since power on
        spec.time_s = spec_timestamp./1000;
        spec.dnum = Meta_Data.start_dnum + days(seconds(spec.time_s));
    end
    
    % Sort epsi fields
    spec = orderfields(spec,{'dnum','freq','time_s','t1_volt','s1_volt','a3_g','data','channels','checksum'});
    catch
        no_data_types = [no_data_types,'spec'];
        spec=[];
    end
end %end spec


%% Process AVGSPEC data
if isempty(ind_avgspec_start)
    no_data_types = [no_data_types,'avgspec'];
    avgspec=[];
else
    processed_data_types = [processed_data_types,'avgspec'];
    %disp('processing avgspec data')
    
    % Pre-allocate space for data
    avgspec.data.n_blocks           = numel(ind_avgspec_start); %Number of data blocks beginning with $EFE
    
    avgspec.data.n_channels        = 3;
    avgspec.data.sample_freq       = 160;
    avgspec.data.sample_per_spec   = 2048/2; %NFFT/2
    avgspec.data.spec_per_block    = 1;
    avgspec.data.timestamp_length  = 8;
    avgspec.data.n_elements        = avgspec.data.timestamp_length+ ...
        avgspec.data.n_channels*...
        avgspec.data.sample_per_spec*4;
    avgspec.channels ={'t1_k','s1_k','a3_g'};
    
    %     seg.data.n_recs             = efe.data.n_blocks*efe.data.recs_per_block;
    % The first 8 bytes are the timestamp. Timestamp in milliseconds since
    % 1970 is too big for uint32, so use uint64
    ind_time = 1:avgspec.data.timestamp_length;
    mult = 256.^[0:7].';
    
    
    avgspec_timestamp           = nan(avgspec.data.n_blocks,1);
    %epsi.time_s and epsi.dnum will be created from epsi_timestamp once
    %all its records are filled
    
    % Loop through data blocks and parse strings
    for iB = 1:avgspec.data.n_blocks
        
        % Grab the block of data starting with the header
        avgspec_block_str = str(ind_avgspec_start(iB):ind_avgspec_stop(iB));
        
        % Get the header timestamp and length of data block (?)
        avgspec.data.hextimestamp.value   = hex2dec(avgspec_block_str(tag.hextimestamp.inds));
        avgspec.data.hexlengthblock.value = hex2dec(avgspec_block_str(tag.hexlengthblock.inds));
        % It should be the case that
        % efe.hexlengthblock.value==efe.data.element_length*efe.data.recs_per_block
        
        % Get the data after the header.
        avgspec_block_data = avgspec_block_str(tag.data_offset:end-tag.chksum.length);
        
        % Does length(efe_block_data)=efe.hexlengthblock.value?
        if (length(avgspec_block_data)~=avgspec.data.n_elements)
            fprintf("AVGSPEC block %i has incorrect length\r\n",iB)
        else
            
            avgspec_timestamp(iB)=double(uint64(avgspec_block_data(ind_time)))*mult;
            
            % Reshape the data
            avgspec.data.raw_bytes{iB} = reshape(uint8(avgspec_block_data(ind_time(end)+1:end)),...
                4,...
                avgspec.data.sample_per_spec*...
                avgspec.data.n_channels);
            
            for n=1:avgspec.data.sample_per_spec*avgspec.data.n_channels
                avgspec.data.raw_float{iB}(n) = double(typecast(avgspec.data.raw_bytes{iB}(:,n).' , 'single'));
            end
            % Save data as 'counts'
            for cha=1:avgspec.data.n_channels
                wh_channel=avgspec.channels{cha};
                avgspec.(wh_channel){iB} = avgspec.data.raw_float{iB}(1+(cha-1)*avgspec.data.sample_per_spec:cha*avgspec.data.sample_per_spec);
            end
            
            % Get the checksum value
            i1 = tag.data_offset+avgspec.data.hexlengthblock.value;
            i2 = i1+tag.chksum.length-1;
            % The full checksum looks like *8F__ so we get str indices from
            % i1+1 to i2-2 to isolate just the 8F part
            avgspec.checksum.data(iB) = hex2dec(avgspec_block_str(i1+1:i2-2));
        end %end if efe data block is the correct size
    end %end loop through efe blocks
    
    [~,fe]=pwelch(1.*(1:Meta_Data.PROCESS.nfft),Meta_Data.PROCESS.nfft,[], ...
        Meta_Data.PROCESS.nfft,Meta_Data.AFE.FS);
    avgspec.freq=fe(2:end);
    
    
    % Clean and concatenate efe.data.raw_bytes
    % Check for empty raw_bytes cells and concatenate only the full ones
    for cha=1:avgspec.data.n_channels
        wh_channel=avgspec.channels{cha};
        idnull=cellfun(@isempty,avgspec.(wh_channel));
        avgspec.(wh_channel)=avgspec.(wh_channel)(~idnull);
    end
    
    if nanmedian(avgspec_timestamp)>1e9
        % time_s - seconds since 1970
        % dnum - matlab datenum
        [avgspec.time_s,avgspec.dnum] = convert_timestamp(avgspec_timestamp);
    else
        % time_s - seconds since power on
        avgspec.time_s = avgspec_timestamp./1000;
        avgspec.dnum = Meta_Data.start_dnum + days(seconds(avgspec.time_s));
        
    end
    
    % Sort epsi fields
    avgspec = orderfields(avgspec,{'dnum','time_s','freq','t1_k','s1_k','a3_g','data','channels','checksum'});
    
end %end spec


%% Process DISRATE data
if isempty(ind_dissrate_start)
    no_data_types = [no_data_types,'dissrate'];
    dissrate=[];
else
    processed_data_types = [processed_data_types,'dissrate'];
    %disp('processing dissrate data')
    
    % Pre-allocate space for data
    dissrate.data.n_blocks           = numel(ind_dissrate_start); %Number of data blocks beginning with $EFE
    
    dissrate.channels         = {'pressure','temperature','salinity','dpdt','chi','chi_fom','epsilon','epsi_fom','nu','kappa'};
    dissrate.data.n_channels         = length(dissrate.channels);
    dissrate.data.dissrate_per_block = 1;
    dissrate.data.timestamp_length   = 8;
    dissrate.data.float_length       = 4;
    dissrate.data.n_elements         = dissrate.data.timestamp_length+ ...
        (dissrate.data.n_channels)*...
        dissrate.data.float_length;
    
    % The first 8 bytes are the timestamp. Timestamp in milliseconds since
    % 1970 is too big for uint32, so use uint64
    ind_time = 1:dissrate.data.timestamp_length;
    mult = 256.^[0:7].';
    
    
    dissrate_timestamp           = nan(dissrate.data.n_blocks,1);
    %epsi.time_s and epsi.dnum will be created from epsi_timestamp once
    %all its records are filled
    
    % Loop through data blocks and parse strings
    for iB = 1:dissrate.data.n_blocks
        
        % Grab the block of data starting with the header
        dissrate_block_str = str(ind_dissrate_start(iB):ind_dissrate_stop(iB));
        
        % Get the header timestamp and length of data block (?)
        dissrate.data.hextimestamp.value   = hex2dec(dissrate_block_str(tag.hextimestamp.inds));
        dissrate.data.hexlengthblock.value = hex2dec(dissrate_block_str(tag.hexlengthblock.inds));
        % It should be the case that
        % efe.hexlengthblock.value==efe.data.element_length*efe.data.recs_per_block
        
        % Get the data after the header.
        dissrate_block_data = dissrate_block_str(tag.data_offset:end-tag.chksum.length);
        
        % Does length(efe_block_data)=efe.hexlengthblock.value?
        if (length(dissrate_block_data)~=dissrate.data.n_elements)
            fprintf("dissrate block %i has incorrect length\r\n",iB)
        else
            
            dissrate_timestamp(iB)=double(uint64(dissrate_block_data(ind_time)))*mult;
            
            % Reshape the data
            dissrate.data.raw_bytes{iB} = reshape(uint8(dissrate_block_data(ind_time(end)+1:end)),...
                dissrate.data.float_length,...
                dissrate.data.n_channels);
            
            for n=1:dissrate.data.n_channels
                dissrate.data.raw_float{iB}(n) = double(typecast(dissrate.data.raw_bytes{iB}(:,n).' , 'single'));
            end
            % Save data as 'counts'
            for cha=1:dissrate.data.n_channels
                wh_channel=dissrate.channels{cha};
                dissrate.(wh_channel){iB} = dissrate.data.raw_float{iB}(cha);
            end
            
            % Get the checksum value
            i1 = tag.data_offset+dissrate.data.hexlengthblock.value;
            i2 = i1+tag.chksum.length-1;
            % The full checksum looks like *8F__ so we get str indices from
            % i1+1 to i2-2 to isolate just the 8F part
            dissrate.checksum.data(iB) = hex2dec(dissrate_block_str(i1+1:i2-2));
        end %end if efe data block is the correct size
    end %end loop through efe blocks
    
    % Clean and concatenate data.raw_bytes
    % Check for empty raw_bytes cells and concatenate only the full ones
    for cha=1:dissrate.data.n_channels
        wh_channel=dissrate.channels{cha};
        idnull=cellfun(@isempty,dissrate.(wh_channel));
        dissrate.(wh_channel)=cell2mat(dissrate.(wh_channel)(~idnull));
    end
    
    if nanmedian(dissrate_timestamp)>1e9
        % time_s - seconds since 1970
        % dnum - matlab datenum
        [dissrate.time_s,dissrate.dnum] = convert_timestamp(dissrate_timestamp);
    else
        % time_s - seconds since power on
        dissrate.time_s = dissrate_timestamp./1000;
        dissrate.dnum = Meta_Data.start_dnum + days(seconds(dissrate.time_s));
        
    end
    
    % Sort epsi fields
    dissrate = orderfields(dissrate,{'dnum','time_s','pressure','temperature', ...
        'salinity','dpdt','chi','chi_fom','epsilon','epsi_fom','nu','kappa',...
        'data','channels','checksum'});
    
end %end spec


%% Process APF data
if (isempty(ind_apf0_start))
    no_data_types = [no_data_types,'apf'];
    apf=[];
else
    processed_data_types = [processed_data_types,'apf'];
    %disp('processing apf data')
    ind_apf_start=ind_apf0_start;
    ind_apf_stop=ind_apf0_stop;
    
    
    sampling_freq=160;
    % Pre-allocate space for data
    apf.data.n_blocks           = numel(ind_apf_start); %Number of data blocks beginning with $EFE
    
    %     apf.data.n_channels         = length(apf.channels);
    apf.data.dissrate_per_block = 1;
    apf.data.avgs_per_block     = 3;
    apf.data.avgs_per_block     = 3;
    apf.data.timestamp_length   = 2;
    apf.data.float_length       = 4;
    apf.data.dissrate_length    = 3;
    apf.data.foco_length        = 2;
    apf.data.fom_length         = 1;
    apf.data.max_dissrate       = log10(1e-3);
    apf.data.max_foco           = log10(1e-2);
    apf.data.max_fom            = 10;
    apf.data.min_dissrate       = log10(1e-12);
    apf.data.min_foco           = log10(1e-13);
    apf.data.min_fom            = 0 ;
    apf.data.dissrate_range     = 4095;
    apf.data.foco_range         = 65535;
    apf.data.fom_range          = 15;
    apf.data.diag_coef          = 8;
    apf.data.dissrate_per_bit   = apf.data.dissrate_range/ ...
        (apf.data.max_dissrate -apf.data.min_dissrate);
    apf.data.foco_per_bit       = apf.data.foco_range/ ...
        (apf.data.max_foco -apf.data.min_foco);
    apf.data.fom_per_bit        = apf.data.fom_range/ ...
        (apf.data.max_fom -apf.data.min_fom);
    
    apf.data.dissrate_count0    = - (apf.data.dissrate_per_bit* ...
        apf.data.max_dissrate) +...
        apf.data.dissrate_range;
    apf.data.foco_count0        = - (apf.data.foco_per_bit* ...
        apf.data.max_foco) +...
        apf.data.foco_range;
    apf.data.fom_count0         = - (apf.data.fom_per_bit* ...
        apf.data.max_fom) +...
        apf.data.fom_range;
    
    %     apf.data.n_elements         = apf.data.timestamp_length+ ...
    %                                       (apf.data.n_channels)*...
    %                                        apf.data.float_length;
    %read the metadata
    %         uint32_t daq_timestamp;                 4//
    %         uint16_t  profile_id;                   2
    %         uint16_t modsom_sn;                     2
    %         uint16_t efe_sn;                        2
    %         uint32_t firmware_rev;                  4
    %         uint16_t nfft;                          2
    %         uint16_t nfftdiag;                      2
    %         mod_som_apf_probe_t  probe1;            5
    %         mod_som_apf_probe_t  probe2;            5
    %         uint8_t  comm_telemetry_packet_format;  1
    %         uint8_t  sd_format;                     1
    %         uint16_t sample_cnt;                    2
    %         uint32_t voltage;                       4
    %         uint16_t end_metadata; //always 0xFFFF; 2
    apf.metadata.size=4+2+2+2+4+2+2+5+5+1+1+2+4+2;
    
    % anonymous function to convert dissrate, foco, fom
    %            mod_epsilon  = (uint32_t) ceil(local_epsilon*
    %            mod_som_apf_ptr->producer_ptr->decim_coef.dissrate_per_bit+
    %            mod_som_apf_ptr->producer_ptr->decim_coef.dissrate_counts_at_origin);
    
    convert_dissrate = @(x) ((x-apf.data.dissrate_count0)/apf.data.dissrate_per_bit);
    convert_foco    = @(x) (x-apf.data.foco_count0)/apf.data.foco_per_bit;
    convert_fom     = @(x) (x-apf.data.fom_count0)/apf.data.fom_per_bit;
    
    
    % The first 8 bytes are the timestamp. Timestamp in milliseconds since
    % 1970 is too big for uint32, so use uint64
    ind_time = 1:apf.data.timestamp_length;
    mult = 256.^[0:7].';
    
    
    %     apf_timestamp           = nan(apf.data.n_blocks,1);
    %epsi.time_s and epsi.dnum will be created from epsi_timestamp once
    %all its records are filled
    
    % Loop through data blocks and parse strings
    for iB = 1:apf.data.n_blocks
        
        % Grab the block of data starting with the header
        apf_block_str = str(ind_apf_start(iB):ind_apf_stop(iB));
        
        % Get the header timestamp and length of data block (?)
        apf.data.hextimestamp.value   = hex2dec(apf_block_str(tag.hextimestamp.inds));
        apf.data.hexlengthblock.value = hex2dec(apf_block_str(tag.hexlengthblock.inds));
        % It should be the case that
        % efe.hexlengthblock.value==efe.data.element_length*efe.data.recs_per_block
        
        % Get the data after the header.
        apf_block_data = apf_block_str(tag.data_offset:end-tag.chksum.length);
        
        %         uint32_t daq_timestamp;                 4//
        %         uint8_t  profile_id;                    1
        %         uint16_t modsom_sn;                     2
        %         uint16_t efe_sn;                        2
        %         uint32_t firmware_rev;                  4
        %         uint16_t nfft;                          2
        %         uint16_t nfftdiag;                      2
        %         mod_som_apf_probe_t  probe1;            5
        %         mod_som_apf_probe_t  probe2;            5
        %         uint8_t  comm_telemetry_packet_format;  1
        %         uint8_t  sd_format;                     1
        %         uint16_t sample_cnt;                    2
        %         uint32_t voltage;                       4
        %         uint16_t end_metadata; //always 0xFFFF; 2
        
        apf.metadata.raw_bytes{iB}     = uint8(apf_block_data(1:apf.metadata.size));
        apf.metadata.daq_timestamp{iB} = double(uint32(apf_block_data(1:4)))*mult(1:4);% 4bytes
        apf.metadata.profile_id{iB}    = double(uint32(apf_block_data(5:6)));% 2 byte
        apf.metadata.modsom_sn{iB}     = double(uint32(apf_block_data(7:8)))*mult(1:2);% 2 bytes
        apf.metadata.efe_sn{iB}        = double(uint32(apf_block_data(9:10)))*mult(1:2);% 2 bytes
        apf.metadata.firmware_rev{iB}  = dec2hex(double(uint32(apf_block_data(11:14)))*mult(1:4));% 4 bytes
        apf.metadata.nfft{iB}          = double(uint32(apf_block_data(15:16)))*mult(1:2);% 2 bytes
        apf.metadata.nfftdiag{iB}      = double(uint32(apf_block_data(17:18)))*mult(1:2);% 2 bytes
        apf.metadata.probe1{iB}.type   = double(uint32(apf_block_data(19)));% 1 bytes
        apf.metadata.probe1{iB}.sn     = double(uint32(apf_block_data(20:21)))*mult(1:2);% 2 bytes
        apf.metadata.probe1{iB}.cal    = double(uint32(apf_block_data(22:23)))*mult(1:2);% 2 bytes
        apf.metadata.probe2{iB}.type   = double(uint32(apf_block_data(24)));% 1 bytes
        apf.metadata.probe2{iB}.sn     = double(uint32(apf_block_data(25:26)))*mult(1:2);% 2 bytes
        apf.metadata.probe2{iB}.cal    = double(uint32(apf_block_data(27:28)))*mult(1:2);% 2 bytes
        apf.metadata.packet_format{iB} = double(uint32(apf_block_data(29)));% 1 byte
        apf.metadata.sd_format{iB}     = double(uint32(apf_block_data(30)));% 1 byte
        apf.metadata.sample_cnt{iB}    = double(uint32(apf_block_data(31:32)))*mult(1:2);% 2 bytes
        apf.metadata.voltage{iB}       = double(uint32(apf_block_data(33:36)))*mult(1:4);% 4 bytes
        apf.metadata.endbytes{iB}      = dec2hex(double(uint32(apf_block_data(37:38)))*mult(1:2));% 2 bytes
        
        
        [~,fe]=pwelch(1.*(1:apf.metadata.nfft{iB}),apf.metadata.nfft{iB},[], ...
            apf.metadata.nfft{iB},sampling_freq);
        apf.data.freq{iB}=fe(2:end);
        apf.data.freq_diag{iB}=fe(1:2:end);
        apf.data.freq_diag{iB}=apf.data.freq_diag{iB}(1:apf.metadata.nfftdiag{iB} );
        
        switch apf.metadata.packet_format{iB}
            case 1
                %time, pressure, dpdt, epsilon, chi, epsi_fom, chi_fom.
                
                apf.channels         = {'dnum', ...
                    'pressure', ...
                    'epsilon','chi','epsi_fom','chi_fom'};
            case 2
                %time, pressure, dpdt, epsilon, chi, avg_t, avg_s, avg_a
                apf.channels         = {'dnum', ...
                    'pressure','temperature','salinity','dpdt', ...
                    'epsilon','chi','kcut_shear','fcut_temp', ...
                    'epsi_fom','chi_fom' ...
                    'avg_tg_k','avg_shear_k','avg_accel_k'};
        end
        if iB==1
            apf.data.n_channels = length(apf.channels);
        end
        
        
        
        % Does length(efe_block_data)=efe.hexlengthblock.value?
        fprintf("apf block length %i\r\n",length(apf_block_data))
        if (length(apf_block_data)~=apf.data.hexlengthblock.value)
            fprintf("apf block %i has incorrect length\r\n",iB)
        else
            
%             local_apf_block_counter=apf.metadata.size+5;
            local_apf_block_counter=apf.metadata.size;
            for d=1:apf.metadata.sample_cnt{iB}
                % get timestamp
                local_apf_block_counter=local_apf_block_counter(end)+ind_time;
                apf.data.timestamp{iB}(d)=(double(apf_block_data(local_apf_block_counter)))*mult(1:2);
                % get pressure
                local_apf_block_counter=local_apf_block_counter(end)+(1:apf.data.float_length);
                apf.data.pressure{iB}(d) = double(typecast(uint8(apf_block_data(local_apf_block_counter)), 'single'));
                if(apf.metadata.packet_format{iB}==2)
                    % get temperature
                    local_apf_block_counter=local_apf_block_counter(end)+(1:apf.data.float_length);
                    apf.data.temperature{iB}(d) = double(typecast(uint8(apf_block_data(local_apf_block_counter)), 'single'));
                    % get salinity
                    local_apf_block_counter=local_apf_block_counter(end)+(1:apf.data.float_length);
                    apf.data.salinity{iB}(d) = double(typecast(uint8(apf_block_data(local_apf_block_counter)), 'single'));
                    % get dpdt
                    local_apf_block_counter=local_apf_block_counter(end)+(1:apf.data.float_length);
                    apf.data.dpdt{iB}(d)     = double(typecast(uint8(apf_block_data(local_apf_block_counter)), 'single'));
                end
                % get dissrate
                local_apf_block_counter=local_apf_block_counter(end)+(1:apf.data.dissrate_length);
                apf.data.raw_dissrate{iB}(d,:) = uint8(apf_block_data(local_apf_block_counter));
                
                tempo_epsilon = bitor( bitshift(uint32(apf.data.raw_dissrate{iB}(d,1)),4), ...
                    bitshift(uint32(apf.data.raw_dissrate{iB}(d,2)),-4));
                apf.data.epsilon{iB}(d) = convert_dissrate(double(tempo_epsilon));
                
                tempo_chi     = bitor( bitshift(uint32(apf.data.raw_dissrate{iB}(d,2)),-8), ...
                    uint32(apf.data.raw_dissrate{iB}(d,3)));
                apf.data.chi{iB}(d) = convert_dissrate(double(tempo_chi));
                
                if(apf.metadata.packet_format{iB}==2)
                    % get kcuttoff shear
                    local_apf_block_counter=local_apf_block_counter(end)+(2:apf.data.float_length+1);
                    apf.data.kcutoff_shear{iB}(d) = double(typecast(uint8(apf_block_data(local_apf_block_counter)), 'single'));
                    % get fcuttoff temp
                    local_apf_block_counter=local_apf_block_counter(end)+(1:apf.data.float_length);
                    apf.data.fcutoff_temp{iB}(d) = double(typecast(uint8(apf_block_data(local_apf_block_counter)), 'single'));
                end
                % get fom
                local_apf_block_counter=local_apf_block_counter(end)+(0:apf.data.fom_length-1);
                apf.data.fom{iB}(d,:) = uint8(apf_block_data(local_apf_block_counter));
                
                tempo_fom_epsi     = uint32(bitshift(apf.data.fom{iB}(d),-4));
                apf.data.epsi_fom{iB}(d) = convert_fom(double(tempo_fom_epsi));
                tempo_fom_chi     = uint32(bitor(apf.data.fom{iB}(d),15));
                apf.data.chi_fom{iB}(d) = convert_fom(double(tempo_fom_chi));
                
                if(apf.metadata.packet_format{iB}==2)
                    for ii=1:apf.metadata.nfftdiag{iB}
                        % get avg sh k spectrum
                        local_apf_block_counter=local_apf_block_counter(end)+(1:apf.data.foco_length);
                        tempo_foco = double(uint8(apf_block_data(local_apf_block_counter)))*mult(1:2);
                        apf.data.avg_shear_k{iB}{d}(ii)=convert_foco(tempo_foco);
                        % get avg tg k spectrum
                        local_apf_block_counter=local_apf_block_counter(end)+(1:apf.data.foco_length);
                        tempo_foco = double(uint8(apf_block_data(local_apf_block_counter)))*mult(1:2);
                        apf.data.avg_tg_k{iB}{d}(ii)=convert_foco(tempo_foco);
                        % get avg a3 k spectrum
                        local_apf_block_counter=local_apf_block_counter(end)+(1:apf.data.foco_length);
                        tempo_foco = double(uint8(apf_block_data(local_apf_block_counter)))*mult(1:2);
                        apf.data.avg_accel_k{iB}{d}(ii)=convert_foco(tempo_foco);
%                         fprintf("d:%i; ii:%i\r\n",d,ii);
                    end
                else
                    % TODO change in APF Firmware the SD write so I do not
                    % have this +1 byte offset
                     local_apf_block_counter=local_apf_block_counter+1;
                end
                
            end %end if efe data block is the correct size
            % Get the checksum value
            i1 = tag.data_offset+apf.data.hexlengthblock.value;
            i2 = i1+tag.chksum.length-1;
            % The full checksum looks like *8F__ so we get str indices from
            % i1+1 to i2-2 to isolate just the 8F part
            apf.checksum.data(iB) = hex2dec(apf_block_str(i1+1:i2-2));
            
            
        end %end loop through efe blocks
    end
    
    % Clean and concatenate data.raw_bytes
    % Check for empty raw_bytes cells and concatenate only the full ones
    for cha=1:apf.data.n_channels
        wh_channel=apf.channels{cha};
        switch wh_channel
            case "dnum"
                dnum=cellfun(@(x,y) datenum(datetime(y+x, 'ConvertFrom', 'posixtime')), ...
                    apf.data.timestamp,apf.metadata.daq_timestamp,'un',0);
                apf.(wh_channel)=cell2mat(dnum);
            case {'pressure','dpdt','fom_epsi','fom_chi'}
                apf.(wh_channel)=cell2mat(apf.data.(wh_channel));
            case {'epsilon','chi'}
                apf.(wh_channel)=10.^cell2mat(apf.data.(wh_channel));
            case {'avg_tg_k','avg_shear_k','avg_accel_k'}
                tempo=cellfun(@(x) cell2mat(x.'),apf.data.(wh_channel),'UniformOutput',false);
                apf.(wh_channel)=10.^cell2mat(tempo.');
        end
    end
    % Sort epsi fields
    switch apf.metadata.packet_format{1}
        case 1
            apf = orderfields(apf,{'channels','metadata','data','checksum','dnum', ...
                'pressure','epsilon','chi'});
        case 2
            apf = orderfields(apf,{'channels','metadata','data','checksum','dnum', ...
                'pressure','dpdt','epsilon','chi', ...
                'avg_tg_k','avg_shear_k','avg_accel_k'});
    end

end %end apf

%% Process APF data
if (isempty(ind_apf2_start))
    no_data_types = [no_data_types,'apf2'];
else
    processed_data_types = [processed_data_types,'apf2'];
    %time, pressure, dpdt, epsilon, chi, avg_t, avg_s, avg_a
    apf.channels         = {'dnum', ...
        'pressure','temperature','salinity','dpdt', ...
        'epsilon','chi','kcut_shear','fcut_temp', ...
        'epsi_fom','chi_fom' ...
        'avg_tg_k','avg_shear_k','avg_accel_k'};
    ind_apf_start=ind_apf2_start;
    ind_apf_stop=ind_apf2_stop;
    
    sampling_freq=160;
    
    % Pre-allocate space for data
    apf.data.n_blocks           = numel(ind_apf_start); %Number of data blocks beginning with $EFE
    
    apf.data.n_channels         = length(apf.channels);
    apf.data.dissrate_per_block = 1;
    apf.data.avgs_per_block     = 3;
    apf.data.avgs_per_block     = 3;
    apf.data.timestamp_length   = 8;
    apf.data.float_length       = 4;
    apf.data.dissrate_length    = [];
    apf.data.foco_length        = [];
    apf.data.fom_length         = [];
    apf.data.max_dissrate       = [];
    apf.data.max_foco           = [];
    apf.data.max_fom            = [];
    apf.data.min_dissrate       = [];
    apf.data.min_foco           = [];
    apf.data.min_fom            = [];
    apf.data.dissrate_range     = [];
    apf.data.foco_range         = [];
    apf.data.fom_range          = [];
    apf.data.diag_coef          = [];
    apf.data.dissrate_per_bit   = [];
    apf.data.foco_per_bit       = [];
    apf.data.fom_per_bit        = [];
    
    apf.data.dissrate_count0    = [];
    apf.data.foco_count0        = [];
    apf.data.fom_count0         = [];
    
    %     apf.data.n_elements         = apf.data.timestamp_length+ ...
    %                                       (apf.data.n_channels)*...
    %                                        apf.data.float_length;
    %read the metadata
    %         uint32_t daq_timestamp;                 4//
    %         uint8_t  profile_id;                    1
    %         uint16_t modsom_sn;                     2
    %         uint16_t efe_sn;                        2
    %         uint32_t firmware_rev;                  4
    %         uint16_t nfft;                          2
    %         uint16_t nfftdiag;                      2
    %         mod_som_apf_probe_t  probe1;            6
    %         mod_som_apf_probe_t  probe2;            6
    %         uint8_t  comm_telemetry_packet_format;  1
    %         uint8_t  sd_format;                     1
    %         uint16_t sample_cnt;                    2
    %         uint32_t voltage;                       4
    %         uint16_t end_metadata; //always 0xFFFF; 2
    apf.metadata.size=4+1+2+2+4+2+2+6+6+1+1+2+4+2;
    
    
    
    % anonymous function to convert dissrate, foco, fom
    %            mod_epsilon  = (uint32_t) ceil(local_epsilon*
    %            mod_som_apf_ptr->producer_ptr->decim_coef.dissrate_per_bit+
    %            mod_som_apf_ptr->producer_ptr->decim_coef.dissrate_counts_at_origin);
    
%     convert_dissrate = @(x) ((x-apf.data.dissrate_count0)/apf.data.dissrate_per_bit);
%     convert_foco    = @(x) (x-apf.data.foco_count0)/apf.data.foco_per_bit;
%     convert_fom     = @(x) (x-apf.data.fom_count0)/apf.data.fom_per_bit;
    
    
    % The first 8 bytes are the timestamp. Timestamp in milliseconds since
    % 1970 is too big for uint32, so use uint64
    ind_time = 1:apf.data.timestamp_length;
    mult = 256.^[0:7].';
    
    
    % ALB Meta data. Do not exist with APF2 yet (01/07/2022)
    %         uint32_t daq_timestamp;                 4//
    %         uint8_t  profile_id;                    1
    %         uint16_t modsom_sn;                     2
    %         uint16_t efe_sn;                        2
    %         uint32_t firmware_rev;                  4
    %         uint16_t nfft;                          2
    %         mod_som_apf_probe_t  probe1;            5
    %         mod_som_apf_probe_t  probe2;            5
    %         uint8_t  comm_telemetry_packet_format;  1
    %         uint8_t  sd_format;                     1
    %         uint16_t sample_cnt;                    2
    %         uint16_t end_metadata; //always 0xFFFF; 2
    %
    %         apf.metadata.raw_bytes{iB}     = uint8(apf_block_data(1:apf.metadata.size));
    %         apf.metadata.daq_timestamp{iB} = double(uint32(apf_block_data(1:4)))*mult(1:4);
    %         apf.metadata.profile_id{iB}    = double(uint32(apf_block_data(5)));
    %         apf.metadata.modsom_sn{iB}     = double(uint32(apf_block_data(6:7)))*mult(1:2);
    %         apf.metadata.efe_sn{iB}        = double(uint32(apf_block_data(8:9)))*mult(1:2);
    %         apf.metadata.firmware_rev{iB}  = dec2hex(double(uint32(apf_block_data(13:16)))*mult(1:4)); % This makes no sense why 13:16 and not 10:13
    %         apf.metadata.nfft{iB}          = double(uint32(apf_block_data(17:18)))*mult(1:2);
    %         apf.metadata.probe1{iB}.type   = double(uint32(apf_block_data(19)));
    %         apf.metadata.probe1{iB}.sn     = double(uint32(apf_block_data(20:21)))*mult(1:2);
    %         apf.metadata.probe1{iB}.cal    = double(uint32(apf_block_data(23:24)))*mult(1:2);% This makes no sense why 23:24 and not 22:23
    %         apf.metadata.probe2{iB}.type   = double(uint32(apf_block_data(25)));
    %         apf.metadata.probe2{iB}.sn     = double(uint32(apf_block_data(26:27)))*mult(1:2);
    %         apf.metadata.probe2{iB}.cal    = double(uint32(apf_block_data(29:30)))*mult(1:2);
    %         apf.metadata.packet_format{iB} = double(uint32(apf_block_data(31)));
    %         apf.metadata.sd_format{iB}     = double(uint32(apf_block_data(32)));
    %         apf.metadata.sample_cnt{iB}    = double(uint32(apf_block_data(33:34)))*mult(1:2);
    %         apf.metadata.endbytes{iB}      = dec2hex(double(uint32(apf_block_data(35:36)))*mult(1:2));
    %         apf.data.nfft_diag             = apf.metadata.nfft{iB}/apf.data.diag_coef;
    
    apf.metadata.nfft=2048;
    [~,fe]=pwelch(1.*(1:apf.metadata.nfft),apf.metadata.nfft,[], ...
        apf.metadata.nfft,sampling_freq);
    
    %epsi.time_s and epsi.dnum will be created from epsi_timestamp once
    %all its records are filled
    
    % Loop through data blocks and parse strings
    for iB = 1:apf.data.n_blocks
        
        % Grab the block of data starting with the header
        apf_block_str = str(ind_apf_start(iB):ind_apf_stop(iB));
        
        % Get the header timestamp and length of data block (?)
        apf.data.hextimestamp.value   = hex2dec(apf_block_str(tag.hextimestamp.inds));
        apf.data.hexlengthblock.value = hex2dec(apf_block_str(tag.hexlengthblock.inds));
        % It should be the case that
        % efe.hexlengthblock.value==efe.data.element_length*efe.data.recs_per_block
        
        % Get the data after the header.
        apf_block_data = apf_block_str(tag.data_offset:end-tag.chksum.length);
        
        apf.data.freq{iB}=fe(2:end);
        
        % Does length(efe_block_data)=efe.hexlengthblock.value?
        if (length(apf_block_data)~=apf.data.hexlengthblock.value && ...
            apf_block_str(apf.data.hexlengthblock.value+33)~='*')
            fprintf("apf block %i has incorrect length\r\n",iB)
        else
            local_apf_block_counter=0;
            % get timestamp
            local_apf_block_counter=local_apf_block_counter(end)+ind_time;
            apf.data.timestamp(iB)=(double(apf_block_data(local_apf_block_counter)))*mult;
            % get pressure
            local_apf_block_counter=local_apf_block_counter(end)+(1:apf.data.float_length);
            apf.data.pressure(iB) = double(typecast(uint8(apf_block_data(local_apf_block_counter)), 'single'));
            % get temperature
            local_apf_block_counter=local_apf_block_counter(end)+(1:apf.data.float_length);
            apf.data.temperature(iB) = double(typecast(uint8(apf_block_data(local_apf_block_counter)), 'single'));
            % get salinity
            local_apf_block_counter=local_apf_block_counter(end)+(1:apf.data.float_length);
            apf.data.salinity(iB) = double(typecast(uint8(apf_block_data(local_apf_block_counter)), 'single'));
            % get dpdt
            local_apf_block_counter=local_apf_block_counter(end)+(1:apf.data.float_length);
            apf.data.dpdt(iB)     = double(typecast(uint8(apf_block_data(local_apf_block_counter)), 'single'));
            % get dissrate
            apf.data.raw_dissrate(iB) = 0;
            local_apf_block_counter=local_apf_block_counter(end)+(1:apf.data.float_length);
            apf.data.epsilon(iB) = double(typecast(uint8(apf_block_data(local_apf_block_counter)), 'single'));
            local_apf_block_counter=local_apf_block_counter(end)+(1:apf.data.float_length);
            apf.data.chi(iB) = double(typecast(uint8(apf_block_data(local_apf_block_counter)), 'single'));
            % get kcuttoff shear
            local_apf_block_counter=local_apf_block_counter(end)+(1:apf.data.float_length);
            apf.data.kcutoff_shear(iB) = double(typecast(uint8(apf_block_data(local_apf_block_counter)), 'single'));
            % get fcuttoff temp
            local_apf_block_counter=local_apf_block_counter(end)+(1:apf.data.float_length);
            apf.data.fcutoff_temp(iB) = double(typecast(uint8(apf_block_data(local_apf_block_counter)), 'single'));
            % get fom
            local_apf_block_counter=local_apf_block_counter(end)+(1:apf.data.float_length);
            apf.data.epsi_fom(iB) = double(typecast(uint8(apf_block_data(local_apf_block_counter)), 'single'));
            local_apf_block_counter=local_apf_block_counter(end)+(1:apf.data.float_length);
            apf.data.chi_fom(iB) = double(typecast(uint8(apf_block_data(local_apf_block_counter)), 'single'));
            
            
            for ii=1:apf.metadata.nfft/2
                % get avg sh k spectrum
                local_apf_block_counter=local_apf_block_counter(end)+(1:apf.data.float_length);
                apf.data.avg_shear_k{iB}(ii) = double(typecast(uint8(apf_block_data(local_apf_block_counter)), 'single'));
            end
            % get avg tg k spectrum
            for ii=1:apf.metadata.nfft/2
                local_apf_block_counter=local_apf_block_counter(end)+(1:apf.data.float_length);
                apf.data.avg_tg_k{iB}(ii) = double(typecast(uint8(apf_block_data(local_apf_block_counter)), 'single'));
            end
            % get avg a3 k spectrum
            for ii=1:apf.metadata.nfft/2
                local_apf_block_counter=local_apf_block_counter(end)+(1:apf.data.float_length);
                apf.data.avg_accel_k{iB}(ii) = double(typecast(uint8(apf_block_data(local_apf_block_counter)), 'single'));
            end
            % Get the checksum value
            i1 = tag.data_offset+apf.data.hexlengthblock.value;
            i2 = i1+tag.chksum.length-1;
            % The full checksum looks like *8F__ so we get str indices from
            % i1+1 to i2-2 to isolate just the 8F part
            apf.checksum.data(iB) = hex2dec(apf_block_str(i1+1:i2-2));
            
            
        end%end loop length(apf_block_data)
    end%end loop through efe blocks
    
    % Clean and concatenate data.raw_bytes
    % Check for empty raw_bytes cells and concatenate only the full ones
    for cha=1:apf.data.n_channels
        wh_channel=apf.channels{cha};
        switch wh_channel
            case "dnum"
                dnum=datenum(datetime(apf.data.timestamp./1000, 'ConvertFrom', 'posixtime'));
                apf.(wh_channel)=dnum;
            case {'pressure','dpdt','temperature','salinity','fom_epsi','fom_chi'}
                apf.(wh_channel)=apf.data.(wh_channel);
            case {'epsilon','chi'}
                apf.(wh_channel)=apf.data.(wh_channel);
            case {'avg_tg_k','avg_shear_k','avg_accel_k'}
                apf.(wh_channel)=cell2mat(apf.data.(wh_channel).');
        end
    end
    % Sort epsi fields
    apf = orderfields(apf,{'channels','metadata','data','checksum','dnum', ...
        'pressure','temperature','salinity','dpdt','epsilon','chi', ...
        'avg_tg_k','avg_shear_k','avg_accel_k'});
end %end apf

%% Process fluoruck data
% $fluor0000000000045cf00000001c*6B0000000000045ceXXXXYYYYZZZZ*57"
% XXXX 4 char bytes to convert to uint_16; 
% YYYY 4 char bytes to convert to uint_16; 
% ZZZZ 4 char bytes to convert to uint_16; 
if isempty(ind_fluor_start)
    no_data_types = [no_data_types,'fluor'];
    fluor=[];
else
    
    sensor.name="tridente";
    sensor.bb.cal0=-1.720756e-6;
    sensor.bb.cal1= 5.871734e4;
    sensor.chla.cal0=-158.44269e-6;
    sensor.chla.cal1= 504.23178e3;
    sensor.fDOM.cal0=-49.728096e-3;
    sensor.fDOM.cal1= 30.125764e3;

%     switch Meta_Data.PROCESS
    sensor.data.sample_freq      = 16;
    sensor.data.sample_period    = 1/sensor.data.sample_freq;
    sensor.data.n_blocks         = numel(ind_fluor_start);
    sensor.data.recs_per_block   = 1;
    sensor.data.n_recs           = sensor.data.n_blocks*sensor.data.recs_per_block;
    sensor.data.timestamp_length = 16;
    sensor.data.length           = 12;   

    processed_data_types = [processed_data_types,'fluor'];

    fluor.data.sample_period    = 1/sensor.data.sample_freq;
    fluor.data.n_blocks         = numel(ind_fluor_start);
    %ALB hard coded number of element per record
    fluor.data.recs_per_block   = 1;
    fluor.data.n_recs           = sensor.data.n_blocks*sensor.data.recs_per_block;
    fluor.data.timestamp_length = 16;
    
    
    % Process fluoruck data
    % --------------------------------
    
    % Pre-allocate space for data
    fluor_timestamp   = nan(sensor.data.n_recs,1);
    %ctd.ctdtime and ctd.ctddnum will be created from ctd_timestamp once
    %all its records are filled
    fluor.bb     = nan(sensor.data.n_recs,1);
    fluor.chla   = nan(sensor.data.n_recs,1);
    fluor.fDOM   = nan(sensor.data.n_recs,1);
    fluor.dnum   = nan(sensor.data.n_recs,1);
    fluor.time_s = nan(sensor.data.n_recs,1);

    % Initialize datarecord counter
    n_rec = 0;

    % Loop through data blocks and parse strings
    for iB = 1:sensor.data.n_blocks

        % Grab the block of data starting with the header
        fluor_block_str = str(ind_fluor_start(iB):ind_fluor_stop(iB));

        % Get the header timestamp and length of data block
        sensor.hextimestamp.value   = hex2dec(fluor_block_str(tag.hextimestamp.inds));
        sensor.hexlengthblock.value = hex2dec(fluor_block_str(tag.hexlengthblock.inds));

        % Get the data after the header.
        fluor_block_data = fluor_block_str(tag.data_offset:end-tag.chksum.length);

        if (length(fluor_block_data)~=sensor.hexlengthblock.value)
            fprintf("fluor block %i has incorrect length\r\n",iB)
        else
            %
            % Where do 16 and 24 come from?
            fluor_block_data=reshape(fluor_block_data,16+sensor.data.length,sensor.data.recs_per_block).';

            for iR=1:sensor.data.recs_per_block

                % Count up the record number
                n_rec=n_rec+1;

                % The ctd data is everything but the last two elements of the
                % block
                element_fluor=fluor_block_data(iR,1:end);

                % The hexadecimal timestamp is the first 16 characters of the
                % data
                fluor_timestamp(n_rec) = hex2dec(element_fluor(1:16));

                % Everything after that is the data
                rec_fluor=element_fluor(sensor.data.timestamp_length+1:end);

                try
                    fluor.bb(n_rec,1)   = hex2dec(rec_fluor(:,1:4));
                    fluor.chla(n_rec,1) = hex2dec(rec_fluor(:,(1:4)+4));
                    fluor.fDOM(n_rec,1) = hex2dec(rec_fluor(:,(1:4)+8));
                catch
                    fluor.bb(n_rec,1)   = nan;
                    fluor.chla(n_rec,1) = nan;
                    fluor.fDOM(n_rec,1) = nan;
                end

                % If timestamp has values like 1.6e12, it is in milliseconds since Jan
                % 1, 1970. Otherwise it's in milliseconds since the start of the record
                if nanmedian(fluor_timestamp)>1e9
                    % time_s - seconds since 1970
                    % dnum - matlab datenum
                    [fluor.time_s(n_rec,1),fluor.dnum(n_rec,1)] = convert_timestamp(fluor_timestamp(n_rec));
                else
                    % time_s - seconds since power on
                    if ~isempty(gps)
                        fluor.dnum(n_rec,1) = gps.dnum(end,1);
                        fluor.time_s(n_rec,1) = (fluor.dnum(n_rec,1)-fluor.dnum(1)).*86400;
                    else
                        fluor.time_s = fluor_timestamp(n_rec)./1000;
                        try % 7/14/24 Tridente time_s is a scalar so this line doesnt work
                        fluor.dnum = Meta_Data.start_dnum + days(seconds(fluor.time_s(n_rec,1)));
                        catch
                            fluor.dnum = Meta_Data.start_dnum;
                        end
                    end
                end

            end
        end %end if fluor data block is the correct size
    end %end loop through fluor blocks
end %end loop if there is fluor data

%% Process ttv data
% $TTVP00000188449d43f100000046*2C00000188449d43f100:28:07 447 ms-000000050 ps,+650 mV,+651 mV,078, 078 *12"
% hh:mm:ss 00:28:07
% aaa        447 ms
% tof       -000000050 ps
% amp_up    +650 mV
% amp_dn    +651 mV
% ratio_up   078
% ratio_down 078

%ALB this is engineering data format.
%ALB will work on other format when they exist
if isempty(ind_ttv_start)
    no_data_types = [no_data_types,'ttv'];
    ttv=[];
else
    
    ttv.data.sample_freq      = 16;
    ttv.data.sample_period    = 1/ttv.data.sample_freq;
    ttv.data.n_blocks         = numel(ind_ttv_start);
    ttv.data.recs_per_block   = 1;
    ttv.data.n_recs           = ttv.data.n_blocks*ttv.data.recs_per_block;
    ttv.data.timestamp_length = 16;
    % ALB I have to change this
    ttv.data.length           = 37;   

    processed_data_types = [processed_data_types,'ttv'];

    
    
    % Process ttv data
    % --------------------------------
    
    % Pre-allocate space for data
    ttv_timestamp   = nan(ttv.data.n_recs,1);
    %ctd.ctdtime and ctd.ctddnum will be created from ctd_timestamp once
    %all its records are filled
    ttv.tof_up   = nan(ttv.data.n_recs,1);
    ttv.tof_down = nan(ttv.data.n_recs,1);
    ttv.dtof     = nan(ttv.data.n_recs,1);
    ttv.vfr      = nan(ttv.data.n_recs,1); %(volume flow rate)
    
    % Initialize datarecord counter
    n_rec = 0;
    
    % Loop through data blocks and parse strings
    for iB = 1:ttv.data.n_blocks

        % Grab the block of data starting with the header
        ttv_block_str = str(ind_ttv_start(iB):ind_ttv_stop(iB));

        % Get the header timestamp and length of data block
        ttv.hextimestamp.value   = hex2dec(ttv_block_str(tag.hextimestamp.inds));
        ttv.hexlengthblock.value = hex2dec(ttv_block_str(tag.hexlengthblock.inds));

        % Get the data after the header.
        ttv_block_data = ttv_block_str(tag.data_offset:end-tag.chksum.length);

        if (length(ttv_block_data)~=ttv.hexlengthblock.value)
            fprintf("TTVP block %i has incorrect length\r\n",iB)
        else
            %
            % Where do 16 and 24 come from?
            ttv_block_data=reshape(ttv_block_data,16+ttv.data.length,ttv.data.recs_per_block).';

            % '00000191536bd52000000000,00000000,00000000,C09DBC00'
            parse_ttv_block_data=sscanf(ttv_block_data.','%016x%08x,%08x,%08x,%08x\r\n');
            if length(parse_ttv_block_data)==ttv.data.recs_per_block*5 % 5 data: timestamp, tofup, tofdown, dtof, vfr (volume flow rate) 

                for iR=1:ttv.data.recs_per_block

                    % Count up the record number
                    n_rec=n_rec+1;

                    % The ctd data is everything but the last two elements of the
                    % block
                    element_ttv=ttv_block_data(iR,1:end);

                    % The hexadecimal timestamp is the first 16 characters of the
                    % data
                    ttv_timestamp(n_rec) = hex2dec(element_ttv(1:16));

                    % Everything after that is the data
                    rec_ttv=element_ttv(ttv.data.timestamp_length+1:end);
                    % parse_ttv_block_data=sscanf(rec_ttv,'%08x,%08x,%08x,%08x\r\n');
                    parse_ttv_block_data=strsplit(rec_ttv,',');

                    % Convert the hexadecimal string to uint32
                    % ttv.tof_up(n_rec)   = uint32(hex2dec(parse_ttv_block_data{1}([7 8 5 6 3 4 1 2])));
                    % ttv.tof_down(n_rec) = uint32(hex2dec(parse_ttv_block_data{2}([7 8 5 6 3 4 1 2])));
                    % ttv.dtof(n_rec)     = uint32(hex2dec(parse_ttv_block_data{3}([7 8 5 6 3 4 1 2])));
                    % ttv.vfr(n_rec)      = uint32(hex2dec(parse_ttv_block_data{4}([7 8 5 6 3 4 1 2])));
                    ttv.tof_up(n_rec)   = double(typecast(uint32(hex2dec(parse_ttv_block_data{1})),'single'));
                    ttv.tof_down(n_rec) = double(typecast(uint32(hex2dec(parse_ttv_block_data{2})),'single'));
                    ttv.dtof(n_rec)     = double(typecast(uint32(hex2dec(parse_ttv_block_data{3})),'single'));
                    ttv.vfr(n_rec)      = double(typecast(uint32(hex2dec(parse_ttv_block_data{4})),'single'));
                    % Convert the uint32 to single precision float
                    % ttv.tof_up(n_rec)   = typecast(ttv.tof_up(n_rec),'single');
                    % ttv.tof_down(n_rec) = typecast(ttv.tof_down(n_rec),'single');
                    % ttv.dtof(n_rec)     = typecast(ttv.dtof(n_rec),'single');
                    % ttv.vfr(n_rec)      = typecast(ttv.tof_vfr(n_rec),'single');


                    % If timestamp has values like 1.6e12, it is in milliseconds since Jan
                    % 1, 1970. Otherwise it's in milliseconds since the start of the record
                    if nanmedian(ttv_timestamp)>1e9
                        % time_s - seconds since 1970
                        % dnum - matlab datenum
                        [ttv.time_s,ttv.dnum] = convert_timestamp(ttv_timestamp);
                    else
                        % time_s - seconds since power on
                        ttv.time_s = ttv_timestamp./1000;
                        ttv.dnum   = Meta_Data.start_dnum + days(seconds(ttv.time_s));
                    end
                end

            end
        end %end if ttv data block is the correct size
    end %end loop through ttv blocks
end %end loop if there is ttv data


fprintf(['   processed data for: ' repmat('%s ',1,length(processed_data_types)), '\n'], processed_data_types{:})
fprintf(['   no data for:        ' repmat('%s ',1,length(no_data_types)), '\n'], no_data_types{:})

% Combine all data
make data epsi ctd alt isap act vnav gps seg spec avgspec dissrate apf fluor ttv Meta_Data


end





%% Timestamp conversion subfunction
% Convert from milliseconds since 1970 to 1) seconds since 1970 and 2)
% Matlab datenum
function [time_s,dnum] = convert_timestamp(timestamp)
seconds1970 = timestamp./1000;
days1970 = seconds1970/(24*60*60);
offset1970_days = datenum(1970,1,1) - datenum(0000,1,0);
offset1970_seconds = offset1970_days*(24*60*60);
time_s = seconds1970 + offset1970_seconds;
dnum = days1970 + offset1970_days;
end

%% SBE calibration subfunctions
function ctd = sbe49_ascii_get_temperature(ctd,sbe)

a0 = sbe.cal.ta0;
a1 = sbe.cal.ta1;
a2 = sbe.cal.ta2;
a3 = sbe.cal.ta3;

mv = (ctd.T_raw-524288)/1.6e7;
r = (mv*2.295e10 + 9.216e8)./(6.144e4-mv*5.3e5);
ctd.T = a0+a1*log(r)+a2*log(r).^2+a3*log(r).^3;
ctd.T = 1./ctd.T - 273.15;
return;
end

%  reads and apply calibration to the conductivity data
function ctd = sbe49_ascii_get_conductivity(ctd,sbe)
try
    g = sbe.cal.g;
    h = sbe.cal.h;
    i = sbe.cal.i;
    j = sbe.cal.j;
    tcor = sbe.cal.tcor;
    pcor = sbe.cal.pcor;
catch
    g = sbe.cal.cg;
    h = sbe.cal.ch;
    i = sbe.cal.ci;
    j = sbe.cal.cj;
    tcor = sbe.cal.ctcor;
    pcor = sbe.cal.cpcor;
end

f = ctd.C_raw/256/1000;

ctd.C = (g+h*f.^2+i*f.^3+j*f.^4)./(1+tcor.*ctd.T+pcor.*ctd.P);

return;
end

%  reads and apply calibration to the pressure data
function ctd = sbe49_ascii_get_pressure(ctd,sbe)
% ALB 04112019 Changed ctd.cal.SBEcal. to ctd.cal.
pa0 = sbe.cal.pa0;
pa1 = sbe.cal.pa1;
pa2 = sbe.cal.pa2;
ptempa0 = sbe.cal.ptempa0;
ptempa1 = sbe.cal.ptempa1;
ptempa2 = sbe.cal.ptempa2;
ptca0 = sbe.cal.ptca0;
ptca1 = sbe.cal.ptca1;
ptca2 = sbe.cal.ptca2;
ptcb0 = sbe.cal.ptcb0;
ptcb1 = sbe.cal.ptcb1;
ptcb2 = sbe.cal.ptcb2;


y = ctd.PT_raw/13107;

t = ptempa0+ptempa1*y+ptempa2*y.^2;
x = ctd.P_raw-ptca0-ptca1*t-ptca2*t.^2;
n = x*ptcb0./(ptcb0+ptcb1*t+ptcb2*t.^2);

ctd.P = (pa0+pa1*n+pa2*n.^2-14.7)*0.689476;

return;
end

%% FCTD header parsing function
%  parse all the lines in the header of the file
function FCTD = FastCTD_ASCII_parseheader(FID)

FCTD = [];
fgetl(FID);
s=fgetl(FID);
[v,val]=FastCTD_ASCII_parseheadline(s);
if ~isempty(v)
    eval(['FCTD.header.' lower(v) '=' val ';']);
end
s=fgetl(FID);
while ~strncmp(s,'%*****END_FCTD',14) && ~feof(FID)
    [v,val]=FastCTD_ASCII_parseheadline(s);
    if ~isempty(v)
        try
            eval(['FCTD.header.' lower(v) '=' val ';']);
        catch obj
            if strncmp(v,'FCTD_VER',8)
                eval(['FCTD.header.' lower(v) '=''' val ''';']);
            else
                if ~contains(v,'$CH')
                    disp(obj.message);
                    disp(['Error occured in string: ' s]);
                end
            end
            
        end
    end
    s=fgetl(FID);
    strncmp(s,'%*****END_FCTD',14);
end
return
end

%  parse each line in the header to detect comments
function [v,val]=FastCTD_ASCII_parseheadline(s)
if s(1)~='%'
    
    i = strfind(s,'=');
    v=s(1:i-1);
    val = s(i+1:end);
else
    v=[];
    val=[];
end

return
end

function cha=parse_epsi_channel_string(str_EPSICHANNEL_header)

    lines_EPSICHANNEL=strsplit(str_EPSICHANNEL_header,'\n');
    cha.ch1=parse_single_epsi_channel(lines_EPSICHANNEL{2});
    cha.ch2=parse_single_epsi_channel(lines_EPSICHANNEL{3});
    cha.ch3=parse_single_epsi_channel(lines_EPSICHANNEL{4});
    cha.ch4=parse_single_epsi_channel(lines_EPSICHANNEL{5});
end

function cha=parse_single_epsi_channel(lines_EPSICHANNEL)
cha.str=lines_EPSICHANNEL;
cha.SN = cha.str(6:8);
calib_str=strsplit(cha.str(11:end),',');
if length(calib_str)==1
    cha.cal=0;
elseif(length(calib_str)==3)
    cha.datecal=calib_str{1};
    cha.cal=str2double(calib_str{2});
    cha.cap=str2double(calib_str{3});
end
end

