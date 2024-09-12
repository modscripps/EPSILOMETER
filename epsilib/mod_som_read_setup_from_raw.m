%read config file
function setup=mod_som_read_setup_from_raw(str)

% fid=fopen(filepath);
% total_str = fread(fid,'*char');
% 
% % NC 13 Aug 2022 - Save filepath so you can open it up again later to get the SBE number
% setup.filepath = filepath;
% 
% [ind_settings_start,ind_settings_stop, ind_settings_tokens] = regexp(total_str.','\$SOM3([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end','tokenExtents');
% 
% % 
% str=total_str(ind_settings_start+32:ind_settings_stop-5);
conv16d=@(x) (x(2).*256+x(1));
conv32d=@(x) (x(4).*256^3+x(3).*256^2+x(2).*256+x(1));
%%
size_offset=1;
size_length=4-1;
%ALB  I am getting the size now so wecan handle previous Firmware  parse config
setup.size=conv32d(double(uint8(str(size_offset+(0:size_length)))));

header_offset=size_offset+size_length+1;
header_length=8-1;
mission_name_offset=header_offset+header_length+1;
mission_name_length=24-1;
vehicle_name_offset=mission_name_offset+mission_name_length+1;
vehicle_name_length=24-1;
firmware_offset=vehicle_name_offset+vehicle_name_length+1;
firmware_length=40-1;

gitid_offset=firmware_offset+firmware_length+1;
switch setup.size
    case {864,896}    
        gitid_length=24-1;
    otherwise
        gitid_length=0-1;
end
rev_offset=gitid_offset+gitid_length+1;
rev_length=8-1;
sn_offset=rev_offset+rev_length+1;
sn_length=8-1;

%ALB I am commenting these lines because
%mod_som_epsi_epsifish2_telemetry-Feb 15 rev3 
% does not have git ID
% TODO change the MODSOM firmware so it has the git ID



% parse config
setup.size=conv32d(double(uint8(str(size_offset+(0:size_length)))));

setup.header=str(header_offset+(0:header_length)).';
setup.header=setup.header(uint8(setup.header)>0);

setup.mission_name=str(mission_name_offset+(0:mission_name_length)).';
setup.mission_name=setup.mission_name(uint8(setup.mission_name)>0);

setup.vehicle_name=str(vehicle_name_offset+(0:vehicle_name_length)).';
setup.vehicle_name=setup.vehicle_name(uint8(setup.vehicle_name)>0);

setup.firmware=str(firmware_offset+(0:firmware_length)).';
setup.firmware=setup.firmware(uint8(setup.firmware)>0);

switch setup.size
    case {864,896}    
        setup.gitid=str(firmware_offset+(0:firmware_length)).';
    otherwise
        setup.gitid=[];
end


initialize_flag_offset=sn_offset+sn_length+1;
initialize_flag_length=4-1;

setup.rev=str(rev_offset+(0:rev_length)).';
setup.rev=setup.rev(uint8(setup.rev)>0);

setup.sn=str(sn_offset+(0:sn_length)).';
setup.sn=setup.sn(uint8(setup.sn)>0);

setup.initialize_flag=conv32d(double(uint8(str(initialize_flag_offset+(0:initialize_flag_length)))));

% TO DO - GET START TIME
setup.start_dnum = 0;

%%
% nb_module=1;
% unpack=true;
% module{nb_module}.index=initialize_flag_offset+4;
% while unpack
%     module{nb_module}.size=conv32d(double(uint8(str(module{nb_module}.index+(0:3)))));
%     module{nb_module}.str=str(module{nb_module}.index+(0:module{nb_module}.size-1));
%     if module{end}.size==0%TDO make it so I am sure length(str) == the end of the modules.
%         unpack=false;
%     else
%         nb_module=nb_module+1;
%         module{nb_module}.index=module{nb_module-1}.index+module{nb_module-1}.size+mod(module{nb_module-1}.size,4);        
%     end
% end


% nb_module=nb_module-1;
%% parse modules
modules_headers=["CALENDAR","CAL", ...
                 "EFE", ...
                 "SBE","SBE49","SBE41","S49","S41","SB49","SB41", ...
                 "SDIO", ...
                 "VOLT", ...
                 "ALTI"];
for i=1:length(modules_headers)
    wh_module=modules_headers{i};
    idx_module=strfind(str,wh_module);
    if ~isempty(idx_module)
        module.size=conv32d(double(uint8(str(idx_module+(-4:-1)))));
        module.str=str(idx_module-4+(0:module.size-1));

        switch wh_module
            case {"CALENDAR","CAL","$CAL"}
                setup.(wh_module)=parse_calendar_module(module);
            case {"SEFE","EFE","EFE3","EFE4"}
                setup.(wh_module)=parse_efe_module(module);
            case "SDIO"
                setup.(wh_module)=parse_sdio_module(module);
            case {"SBE49","SBE","S49","SB49"}
                setup.(wh_module)=parse_sbe49_module(module);
            case {"SBE41","S41","SB41"}
                setup.(wh_module)=parse_sbe49_module(module);
            case {"ALT","ALTI"}
                setup.(wh_module)=parse_altimeter_module(module);
            case {"VOLT","VOL"}
                setup.(wh_module)=parse_voltage_module(module);
        end
    end
end

% for i=1:nb_module
%     id_module=cellfun(@(x)(strfind(module{i}.str,x)),modules_headers,'un',0);
%     id_module=~cellfun(@isempty,id_module);
%     try
%         wh_module=modules_headers{id_module};
%         setup.(wh_module).str=module{i}.str;
%         switch wh_module
%             case {"CALENDAR","CAL","$CAL"}
%                 setup.(wh_module)=parse_calendar_module(setup.(wh_module));
%             case {"SEFE","EFE","EFE3","EFE4"}
%                 setup.(wh_module)=parse_efe_module(setup.(wh_module));
%             case "SDIO"
%                 setup.(wh_module)=parse_sdio_module(setup.(wh_module));
%             case {"SBE49","SBE","S49","SB49"}
%                 setup.(wh_module)=parse_sbe49_module(setup.(wh_module));
%             case {"SBE41","S41","SB41"}
%                 setup.(wh_module)=parse_sbe49_module(setup.(wh_module));
%             case {"ALT","ALTI"}
%                 setup.(wh_module)=parse_altimeter_module(setup.(wh_module));
%             case {"VOLT","VOL"}
%                 setup.(wh_module)=parse_voltage_module(setup.(wh_module));
%         end
%     catch
%         fprintf("No module name for %s.\r",module{i}.str(1:10))
%     end
% end %end loop through modules
% 
end %end mod_som_read_setup_from_raw

%% CALENDAR
% typedef struct{
% 	uint32_t size;
% 	char header[MOD_SOM_CALENDAR_MAXIMUM_LENGTH_HEADER]; //maximum
%     sl_sleeptimer_date_t initial_date; //Watch for the settings offsets. the size here is not a multiple of words.
%     uint32_t initialize_flag;
% }mod_som_calendar_settings_t,*mod_som_calendar_settings_ptr_t;
%
% typedef  struct  time_date {
%   uint8_t sec;
%   uint8_t min;
%   uint8_t hour;
%   uint8_t month_day;
%   sl_sleeptimer_month_t month;    //uint8_t
%   uint16_t year;                              ///< Current year, based on a 1900 Epoch.
%   sl_sleeptimer_weekDay_t day_of_week; //uint8_t
%   uint16_t day_of_year;
%   sl_sleeptimer_time_zone_offset_t time_zone; ///< Offset, in seconds,
%   from UTC //uint32_t
% } sl_sleeptimer_date_t;

% function calendar=parse_calendar_str(calendar)
function calendar=parse_calendar_module(calendar)

% conv16d=@(x) (x(2).*256+x(1));
conv32d=@(x) (x(4).*256^3+x(3).*256^2+x(2).*256+x(1));
size_offset=1;
size_length=4-1;
header_offset=size_offset+size_length+1;
header_length=16-1;
date_offset=header_offset+header_length+1;
date_length=15;
initialized_flag_offset=date_offset+date_length+1;
initialized_flag_length=4-1;
% parse config
calendar.size=conv32d(double(uint8(calendar.str(size_offset+(0:size_length)))));
calendar.header=calendar.str(header_offset+(1:header_length)).';
calendar.date=parse_date_str(calendar.str(date_offset+(0:date_length)));
calendar.initialized_flag=conv32d(double(uint8(calendar.str(initialized_flag_offset+(0:initialized_flag_length)))));

end

%% voltage module
% *******************************************************************************
%  * @brief
%  *   voltage module setting
%  *
%  *
%  ******************************************************************************/
% typedef struct{
%   uint32_t size;
%   char header[8];
%   uint32_t sampling_frequency;
%   uint32_t initialize_flag;
% 
% }mod_som_voltage_settings_t,*mod_som_voltage_settings_ptr_t;

function volt=parse_voltage_module(volt)

conv16d=@(x) (x(2).*256+x(1));
conv32d=@(x) (x(4).*256^3+x(3).*256^2+x(2).*256+x(1));
size_offset=1;
size_length=4-1;
header_offset=size_offset+size_length+1;
header_length=8-1;
sampling_offset=header_offset+header_length+1;
sampling_length=4-1;
initialized_flag_offset=sampling_offset+sampling_length+1;
% parse config
volt.size=conv32d(double(uint8(volt.str(size_offset+(0:size_length)))));
volt.header=volt.str(header_offset+(1:header_length)).';
volt.sampling=conv32d(double(volt.str(sampling_offset+(0:sampling_length))));
volt.initialized_flag=uint8(volt.str(initialized_flag_offset));

end

%% altimeter module
% typedef struct{
%     uint16_t size;
%     char header[3];
%     char rev[5];
%     char sn[3];
%     uint32_t tx_repetition_period; // see details above
%     uint8_t  tx_repetition_mode;   // see details above
%     uint32_t tx_pulse_width;       // micro-seconds
%     uint32_t blanking_interval;    // micro-seconds
%     bool initialize_flag;
% }
% mod_som_altimeter_settings_t, *mod_som_altimeter_settings_ptr_t;


function altimeter=parse_altimeter_module(altimeter)

conv16d=@(x) (x(2).*256+x(1));
conv32d=@(x) (x(4).*256^3+x(3).*256^2+x(2).*256+x(1));
size_offset=1;
size_length=2-1;
header_offset=size_offset+size_length+1;
header_length=3-1;
rev_offset=header_offset+header_length+1;
rev_length=5-1;
sn_offset=rev_offset+rev_length+1;
sn_length=3-1;
tx_repetition_period_offset=sn_offset+sn_length+1+3;
tx_repetition_period_length=4-1;
tx_repetition_mode_offset=tx_repetition_period_offset+tx_repetition_period_length+1;
tx_repetition_mode_length=1-1;
tx_pulse_width_offset=tx_repetition_mode_offset+tx_repetition_mode_length+1+3;
tx_pulse_width_length=4-1;
blanking_interval_offset=tx_pulse_width_offset+tx_pulse_width_length+1;
blanking_interval_length=4-1;

% parse config
altimeter.size=conv16d(double(uint8(altimeter.str(size_offset+(0:size_length)))));
altimeter.header=altimeter.str(header_offset+(0:header_length)).';
altimeter.rev=altimeter.str(rev_offset+(0:rev_length)).';
altimeter.sn=altimeter.str(sn_offset+(0:sn_length)).';
altimeter.tx_repetition_period=conv32d(double(altimeter.str(tx_repetition_period_offset+(0:tx_repetition_period_length))));
altimeter.tx_repetition_mode=double(altimeter.str(tx_repetition_mode_offset+(0:tx_repetition_mode_length)));
altimeter.tx_pulse_width=conv32d(double(altimeter.str(tx_pulse_width_offset+(0:tx_pulse_width_length))));
altimeter.blanking_interval=conv32d(double(altimeter.str(blanking_interval_offset+(0:blanking_interval_length))));
end



%% EFE
% typedef struct{
%     uint32_t size;
%     char header[8];
%     char rev[8];
%     char sn[8];
%     uint32_t number_of_channels;
%     uint32_t nb_sample_per_record;
%     uint32_t nb_record_per_buffer;
%     uint32_t spi_baudrate;
%     sensor_spec_t sensors[7];
%     uint32_t initialize_flag;
% }
% mod_som_efe_setup_t, *mod_som_efe_setup_ptr_t;

% typedef struct sensor_spec_t{
% 	char			 	name[4]; 		// Sensor Name/ID
% 	char			 	sn[4]; 		    // Sensor serial Number
% 	float		    	cal; 		    // cal is the calibration coefficient
% 	AD7124 				registers; 		// Software register emulators
% 	uint32_t			selector_cs_code; 			// Location
% } sensor_spec_t, *sensor_spec_ptr_t;

% typedef struct AD7124_REGISTERS {
% 	uint8_t 	COMMS;
% 	uint8_t		STATUS;
% 	uint16_t	ADC_CONTROL;
% 	uint32_t 	DATA;
% 	uint32_t 	IO_CONTROL_1;
% 	uint16_t	IO_CONTROL_2;
% 	uint8_t		ID;
% 	uint32_t	ERROR;
% 	uint32_t	ERROR_EN;
% 	uint8_t		MCLK_COUNT;
% 	uint16_t	CHANNEL_0;
% 	uint16_t	CONFIG_0;
% 	uint32_t	FILTER_0;
% 	uint32_t	OFFSET_0;
% 	uint32_t	GAIN_0;
% } AD7124;

function efe=parse_efe_module(efe)
conv16d=@(x) (x(2).*256+x(1));
conv32d=@(x) (x(4).*256^3+x(3).*256^2+x(2).*256+x(1));

size_offset   = 1;
size_length   = 4-1;
header_offset = size_offset+size_length+1;
header_length = 8-1; %it is "$EFE" currently TODO change the firmware to EFE 
rev_offset    = header_offset+header_length+1;
rev_length    = 8-1;
sn_offset     = rev_offset+rev_length+1;
sn_length     = 8-1;

nb_channel_offset=sn_offset+sn_length+1; % WARNING I do not understand why I have to add 2 instead of 1. Maybe a padding issue
nb_channel_length=4-1;
sample_per_record_offset=nb_channel_offset+nb_channel_length+1;
sample_per_record_length=4-1;
nb_record_per_buffer_offset=sample_per_record_offset+sample_per_record_length+1;
nb_record_per_buffer_length=4-1;
spi_baudrate_offset=nb_record_per_buffer_offset+nb_record_per_buffer_length+1;
spi_baudrate_length=4-1;
sensor_offset=spi_baudrate_offset+spi_baudrate_length+1;
sensor_length=60-1;

efe.size=conv32d(double(uint8(efe.str(size_offset+(0:size_length)))));
efe.header=efe.str(header_offset+(0:header_length)).';
efe.header=efe.header(uint8(efe.header)>0);

efe.rev=efe.str(rev_offset+(0:rev_length)).';
efe.rev=efe.rev(uint8(efe.rev)>0);
efe.sn=efe.str(sn_offset+(0:sn_length)).';
efe.sn=efe.sn(uint8(efe.sn)>0);

efe.nb_channel=conv32d(double(efe.str(nb_channel_offset+(0:nb_channel_length))));
efe.sample_per_record=conv32d(double(efe.str(sample_per_record_offset+(0:sample_per_record_length))));
efe.nb_record_per_buffer=conv32d(double(efe.str(nb_record_per_buffer_offset+(0:nb_record_per_buffer_length))));
efe.spi_baudrate=conv32d(double(efe.str(spi_baudrate_offset+(0:spi_baudrate_length))));

for i=1:efe.nb_channel
    efe.sensors{i}=parse_efe_sensor(efe.str( ( ((i-1)*(sensor_length+1))+sensor_offset+(0:sensor_length))));
end

initialized_flag_offset=sensor_offset+efe.nb_channel*(sensor_length+1);
efe.initialized_flag=uint8(efe.str(initialized_flag_offset));
end

%% SBE49
% typedef struct{
%     uint32_t size;
%     char data_header_text[8];
%     char sn[8];
%     uint32_t data_format;
%     uint32_t elements_per_record; //i.e. stream/store frequency.
%     uint32_t initialize_flag;
% }mod_som_sbe49_setup_t, *mod_som_sbe49_setup_ptr_t;

function sbe49=parse_sbe49_module(sbe49)
conv16d=@(x) (x(2).*256+x(1));
conv32d=@(x) (x(4).*256^3+x(3).*256^2+x(2).*256+x(1));

MOD_SOM_SBE49_OUTPUT0_SAMPLE_DATA_LENGTH=24;
MOD_SOM_SBE49_OUTPUT1_SAMPLE_DATA_LENGTH=24;
MOD_SOM_SBE49_OUTPUT2_SAMPLE_DATA_LENGTH=24;
MOD_SOM_SBE49_OUTPUT3_SAMPLE_DATA_LENGTH=40;


size_offset=1;
size_length=4-1;
header_offset=size_offset+size_length+1;
header_length=8-1; %it is "EFE" currently TODO change the firmware to EFE 
sn_offset=header_offset+header_length+1;
sn_length=8-1;
data_format_offset=sn_offset+sn_length+1;
data_format_length=4-1;
sample_data_per_record_offset=data_format_offset+data_format_length+1;
sample_data_per_record_length=4-1;
initialize_flag_offset=sample_data_per_record_offset+sample_data_per_record_length+1;
initialize_flag_length=4-1;

sbe49.size=conv32d(double(uint8(sbe49.str(size_offset+(0:size_length)))));
sbe49.header=sbe49.str(header_offset+(0:header_length)).';
sbe49.header=sbe49.header(uint8(sbe49.header)>0);
sbe49.sn=sbe49.str(sn_offset+(0:sn_length)).';
sbe49.sn=sbe49.sn(uint8(sbe49.header)>0);
sbe49.data_format=conv32d(double(uint8(sbe49.str(data_format_offset+(0:data_format_length)))));
switch sbe49.data_format
    case 0
        sbe49.sample_data_length=MOD_SOM_SBE49_OUTPUT0_SAMPLE_DATA_LENGTH;
    case 1
        sbe49.sample_data_length=MOD_SOM_SBE49_OUTPUT1_SAMPLE_DATA_LENGTH;
    case 2
        sbe49.sample_data_length=MOD_SOM_SBE49_OUTPUT2_SAMPLE_DATA_LENGTH;
    case 3
        sbe49.sample_data_length=MOD_SOM_SBE49_OUTPUT3_SAMPLE_DATA_LENGTH;
end
sbe49.sample_data_per_record=conv32d(double(sbe49.str(sample_data_per_record_offset+(0:sample_data_per_record_length))));
sbe49.initialize_flag=conv32d(double(uint8(sbe49.str(initialize_flag_offset+(0:initialize_flag_length)))));
end



%% SDIO
% typedef struct{
% 	uint32_t size;
% 	char  header[8];
% 	uint32_t file_duration;
% 	char prefix_file[40];
% 	uint32_t initialize_flag;
% }mod_som_sdio_setup_t, *mod_som_sdio_setup_ptr_t;

function sdio=parse_sdio_module(sdio)
conv16d=@(x) (x(2).*256+x(1));
conv32d=@(x) (x(4).*256^3+x(3).*256^2+x(2).*256+x(1));

size_offset=1;
size_length=4-1;
header_offset          = size_offset+size_length+1;
header_length          = 8-1; %it is "$EFE" currently TODO change the firmware to EFE 
file_duration_offset   = header_offset+header_length+1;
file_duration_length   = 4-1;
prefix_file_offset     = file_duration_offset+file_duration_length+1;
prefix_file_length     = 40-1;
initialize_flag_offset = prefix_file_offset+prefix_file_length+1;
initialize_flag_length = 4-1;

sdio.size=conv32d(double(uint8(sdio.str(size_offset+(0:size_length)))));
sdio.header=sdio.str(header_offset+(0:header_length)).';
sdio.header=sdio.header(uint8(sdio.header)>0);
sdio.file_duration=conv32d(double(uint8(sdio.str(file_duration_offset+(0:file_duration_length)))));
sdio.prefix_file=sdio.str(prefix_file_offset+(0:prefix_file_length)).';
sdio.prefix_file=sdio.prefix_file(sdio.prefix_file~=sdio.prefix_file(end));
sdio.initialize_flag=conv32d(double(uint8(sdio.str(initialize_flag_offset+(0:initialize_flag_length)))));
end


%% functions 
function date_matlab=parse_date_str(str)
offset=1;
conv16d=@(x) (x(2).*256+x(1));
conv32d=@(x) (x(4).*256^3+x(2).*256^2+x(2).*256+x(1));
orig_year=1900;

date_t.sec=uint8(str(offset));
date_t.min=uint8(str(offset+1));
date_t.hour=uint8(str(offset+2));
date_t.monthday=uint8(str(offset+3));
date_t.month=uint8(str(offset+4));
date_t.year=orig_year+conv16d(double(uint8(str(offset+(6:7)))));
date_t.dayofweek=uint8(str(offset+8));
date_t.dayofyear=conv16d(double(uint8(str(offset+(10:11)))));
date_t.time_zone=conv32d(double(uint8(str(offset+(11:15)))));

date_str=sprintf("%04i/%02i/%02i %02i:%02i:%02i", ...
                 date_t.year,date_t.month,date_t.monthday+1,...
                 date_t.hour,date_t.min,date_t.sec);
date_matlab=datenum(date_str);             
end


function sensor=parse_efe_sensor(str)
conv16d=@(x) (x(2).*256+x(1));
conv32d=@(x) (x(4).*256^3+x(3).*256^2+x(2).*256+x(1));

sensor_names={'t1','t2','s1','s2','a1','a2','a3','fluor','cond'};
% sensor.name = str(offset+(0:2)).';
% sensor.name = sensor.name(uint8(sensor.name)>0);%str can be paded with null (uint8 0);
for n=1:length(sensor_names)
    wh_name=sensor_names{n};
    idx_name=strfind(str,wh_name);
    if~isempty(idx_name)
        offset=idx_name;
        sensor.name=str(idx_name+(0:3));
        sensor.name=sensor.name(uint8(sensor.name)>0);
    end
end
sensor.sn   = str(offset+(4:6)).';
offset=offset+2;
sensor.cal   = typecast( uint8(str(offset+(6:9))) , 'single');
sensor.register.COMMS         = dec2hex(uint8(str(offset+(10))));
sensor.register.STATUS        = dec2hex(uint8(str(offset+(11))));
% offset=offset+2;
sensor.register.ADC_CONTROL   = dec2hex(conv16d(double(uint8(str(offset+(12:13))))));
sensor.register.DATA          = dec2hex(conv32d(double(uint8(str(offset+(14:17))))));
sensor.register.IO_CONTROL1   = dec2hex(conv32d(double(uint8(str(offset+(18:21))))));
sensor.register.IO_CONTROL2   = dec2hex(conv16d(double(uint8(str(offset+(22:23))))));
sensor.register.ID            = dec2hex(uint8(str(offset+(24))));
sensor.register.ERROR         = dec2hex(conv32d(double(uint8(str(offset+(25:28))))));
sensor.register.ERROR_EN      = dec2hex(conv32d(double(uint8(str(offset+(29:32))))));
sensor.register.MCLK_COUNT    = dec2hex(uint8(str(offset+(33))));
offset=offset+2;
sensor.register.CHANNEL_0     = dec2hex(conv16d(double(uint8(str(offset+(34:35))))));
sensor.register.CONFIG_0      = dec2hex(conv16d(double(uint8(str(offset+(36:37))))));
offset=offset+2;
sensor.register.FILTER_0      = dec2hex(conv32d(double(uint8(str(offset+(38:41))))));
sensor.register.OFFSET_0      = dec2hex(conv32d(double(uint8(str(offset+(42:45))))));
sensor.register.GAIN_0        = dec2hex(conv32d(double(uint8(str(offset+(46:49))))));
sensor.selector_cs_code       = dec2hex(conv32d(double(uint8(str(offset+(50:53))))));
end