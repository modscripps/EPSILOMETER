function epsi = mod_read_epsi_raw(filename,Meta_Data)
% SN_READ_EPSI_RAW - reads epsi raw data
%
% SN_READ_EPSI_RAW(FILENAME) returns a EPSI structure of variables described
% for EPSI data files
% 
% SN_READ_EPSI_RAW(DIRNAME) reads all *.epsi, *.dat, *.bin files in the directory DIRNAME
% 
% SN_READ_EPSI_RAW({FILE1, FILE2,...}) reads all files indicated
%
% this function does not sort data yet nor does it have data checksum
%
% Created 2018/10/15 San Nguyen
%
% add Meta_Data to automatize unipolar and bipolar config of the ADC and the
% number channels
% modified 2018/12/20 Arnaud Le Boyer 
%
% 10/15/20 modifications by Nicole Couto:
%   - rename t1, s1, a1, etc to t1_volts, s1_volts, a1_g

% Check if it is a single file or a directory and a set of files. If it's
% one file, open it. If it's several, combine them into one file. Send
% whatever you have back through this same function until you have an open
% file that will be read by the subfunction mod_read_epsi_raw_file
if ischar(filename) % dir or file
    switch exist(filename,'file')
        case 2 % if it is a file
            fid = fopen(filename,'r');
            epsi = mod_read_epsi_raw(fid,Meta_Data);
            fclose(fid);
        case 7 % if it is a directory
            my_epsi_file = [];
            my_epsi_file = [my_epsi_file; dir(fullfile(filename, '*.epsi'))];
            my_epsi_file = [my_epsi_file; dir(fullfile(filename, '*.bin'))];
            my_epsi_file = [my_epsi_file; dir(fullfile(filename, '*.dat'))];
            if isempty(my_epsi_file)
                epsi = [];
                return
            else
                % prepare to read all files
                epsi = cell(size(my_epsi_file));
                % read the files in the directory
                for i = 1:length(my_epsi_file)
                    disp(['reading ' my_epsi_file(i).name]);
                    epsi{i} = mod_read_epsi_raw(fullfile(filename,my_epsi_file(i).name),Meta_Data);
                end
                % combine all files into one MET structure
                epsi = mod_combine_epsi(epsi{:});
            end
        otherwise
            error('MATLAB:mod_read_epsi_raw:wrongFILENAME','Invalid file specified.');
    end
elseif iscellstr(filename) % cell of files
    % prepare to read all files
    epsi = cell(size(filename));
    % read all files
    for i = 1:length(filename)
        disp(['reading ' filename{i}]);
        epsi{i} = mod_read_epsi_raw(filename{i},Meta_Data);
    end
    % combine all files into one epsi structure
    epsi = mod_combine_epsi(epsi{:});
else
    
    if (filename<1)
        error('MATLAB:mod_read_epsi_raw:wrongFID','FID is invalid');
    end
    
    epsi = mod_read_epsi_raw_file(filename,Meta_Data);
    
    return
end

end

% reading epsi files through FID
function EPSI = mod_read_epsi_raw_file(fid,Meta_Data)
%tic
% make sure the file marker begins at the start
EPSI = epsi_ascii_parseheader(fid);

if ~isempty(EPSI)
    %convert time to MATLAB time
    if EPSI.header.offset_time < 0
        EPSI.header.offset_time = epsi_ascii_correct_negative_time(EPSI.header.offset_time)/86400+datenum(1970,1,1);
    else
        EPSI.header.offset_time = EPSI.header.offset_time/86400+datenum(1970,1,1);
    end
    if EPSI.header.system_time < 0
        EPSI.header.system_time = epsi_ascii_correct_negative_time(EPSI.header.system_time)/86400/100+EPSI.header.offset_time;
    else
        EPSI.header.system_time = EPSI.header.system_time/86400/100+EPSI.header.offset_time;
    end
else
%  isfield(Meta_Data,'SBEcal')
%     EPSI.header=Meta_Data.SBEcal;
    EPSI.header=Meta_Data.aux1.cal;
end

%TODO get rid of ths offset  
switch Meta_Data.PROCESS.recording_mode
    case 'STREAMING'
    offset1=13;
    case 'SD'
    offset1=2;
end


fseek(fid,0,1);
fsize = ftell(fid);
frewind(fid);

str = fread(fid,'*char')';



%
%clc;
% tic
% get MADRE position beginning of datablock
ind_madre = strfind(str,'$MADRE');
% get aux1 position beginning of SBE block if present
%ind_aux1 = strfind(str,'$AUX1');
is_aux1 = contains(str,'$AUX1');
%toc

%find header length from the 1st MADRE block
header_length=strfind(str(ind_madre(1):ind_madre(1)+72),'$');
header_length=header_length(2);
switch header_length
    case 71 % ALB: I added the microsecond timestamp in third position. So the header is longer
        firmware_version='microsecond';
    otherwise
        firmware_version='SODA';
end


% hard coded value to define MADRE block Header
madre.epsi_stamp_length = 8;
madre.epsi_time_length = 8;
madre.alt_time_length = 4;
madre.aux_chksum_length = 8;
madre.map_chksum_length = 8;
madre.fsync_err_length = 8;
madre.offset = 0;
madre.name_length = 6;

switch firmware_version
    case 'microsecond' % ALB: I added the microsecond timestamp in third position. So the header is longer
        madre.epsi_mutime_length = 8;
        madre.epsi_stamp_offset  = -1+madre.name_length;
        madre.epsi_time_offset   = madre.epsi_stamp_offset+madre.epsi_stamp_length+1;
        madre.epsi_mutime_offset = madre.epsi_time_offset+madre.epsi_time_length+1;
        madre.fsync_err_offset    = madre.epsi_mutime_offset+madre.epsi_time_length+1;
        madre.aux_chksum_offset  = madre.fsync_err_offset+madre.fsync_err_length+1;
        madre.alt_time_offset    = madre.aux_chksum_offset+madre.aux_chksum_length+[0 4]+1;
        madre.map_chksum_offset  = madre.alt_time_offset(end)+madre.alt_time_length+1;
    otherwise
        % hard coded value to define MADRE block Header
       
        madre.epsi_stamp_offset  = -1+madre.name_length;
        madre.epsi_time_offset   = madre.epsi_stamp_offset+madre.epsi_stamp_length+1;
        madre.fsync_err_offset    = madre.epsi_time_offset+madre.epsi_time_length+1;
        madre.aux_chksum_offset  = madre.fsync_err_offset+madre.fsync_err_length+1;
        madre.alt_time_offset    = madre.aux_chksum_offset+madre.aux_chksum_length+[0 4]+1;
        madre.map_chksum_offset  = madre.alt_time_offset(end)+madre.alt_time_length+1;
        
end

madre.header_offset  = madre.map_chksum_offset+madre.map_chksum_length+2;

% define offset if aux1 is present
if is_aux1
    if strcmp(Meta_Data.mission,"PILOT2018") && ...
       (strcmp(Meta_Data.deployment,"d2")    || ...
       strcmp(Meta_Data.deployment,"d3")    || ...
       strcmp(Meta_Data.deployment,"d4"))
            aux1.name_length  = 5+2; % +2 becasue there is \r\n for inside the drop files
            aux1.stamp_length = 8; % length of epsi sample number linked to SBE sample.
            aux1.sbe_length   = 22;  % length of SBE sample.
            aux1.nbsample     = 9;
    else
            aux1.name_length  = 5;
            aux1.stamp_length = 8; % length of epsi sample number linked to SBE sample.
            aux1.sbe_length   = 22;  % length of SBE sample.
            aux1.nbsample     = 9;
    end
    aux1_sample_length= aux1.stamp_length + 1 + aux1.sbe_length + 2;
    % e.g 00000F2E,052C2409E6F3080D7A4DAF
    aux1.stamp_offset  = madre.header_offset + (0:aux1.nbsample-1)*aux1_sample_length+aux1.name_length;
    aux1.sbe_offset    = madre.header_offset + (0:aux1.nbsample-1)*aux1_sample_length+(aux1.stamp_length + 1)+aux1.name_length;
end

if is_aux1
    epsi.offset = aux1.stamp_offset+aux1_sample_length;
else
    epsi.offset = madre.header_offset;
end



epsi.name_length = 5; % length of epsi block header ($EPSI).
epsi.nbsamples = 160;   % number of sample in 1 epsi block.
epsi.nchannels = Meta_Data.PROCESS.nb_channels; % number of channels defined by user.
epsi.sample_freq = 320; % hardcoded sampling rate can do better when we will use usecond resolution for the timer. 
epsi.sample_period = 1/epsi.sample_freq;
epsi.bytes_per_adc=3; % data in hex
epsi.bytes_per_channel=epsi.bytes_per_adc; % data in hex
epsi.bytes_per_channel = 3; % length of an ADC sample = length of one channel sample. In bytes MSBF (TO DO check MSBF) 

hexepsi=0;
nb_bytes_end_sample=0;
if strcmp(Meta_Data.mission,"PILOT2018") && ...
   (strcmp(Meta_Data.deployment,"d2")    || ...
   strcmp(Meta_Data.deployment,"d3")    || ...
   strcmp(Meta_Data.deployment,"d4"))
    epsi.offset=epsi.offset+2;
    hexepsi=1;
    epsi.bytes_per_channel=epsi.bytes_per_adc*2; % data in hex
    nb_bytes_end_sample=2;
end
epsi.total_length = epsi.nbsamples*(epsi.nchannels*epsi.bytes_per_channel+nb_bytes_end_sample); % length of an EPSI block




% find the non corrupted (right length)
%!!!!!!!! VERY IMPORTANT to remember !!!!! 
indblock=ind_madre(diff(ind_madre)==median(diff(ind_madre))); %TO DO figure that exact math for the block length
%indblock=ind_madre(diff(ind_madre)==epsi.offset(end)+epsi.name_length+epsi.total_length+offset1+1);
NBblock=numel(indblock);

% initialize arrays and structures.
if(isfield(EPSI,'header'))
    system.time = char(zeros(NBblock,11));
end

madre.epsi_stamp = char(zeros(NBblock,madre.epsi_stamp_length));
madre.epsi_time = char(zeros(NBblock,madre.epsi_time_length));
madre.altimeter = char(zeros(NBblock*2,madre.alt_time_length));
madre.fsync_err = char(zeros(NBblock,madre.fsync_err_length));
madre.aux1_chksum = char(zeros(NBblock,madre.aux_chksum_length));
madre.epsi_chksum = char(zeros(NBblock,madre.map_chksum_length));

switch firmware_version
    case 'microsecond' % ALB: I added the musecond timestamp in third position. So the header is longer
        madre.epsi_mutime = char(zeros(NBblock,madre.epsi_time_length));
        EPSI.madre = struct(...
            'EpsiStamp',NaN(NBblock,1),...
            'TimeStamp',NaN(NBblock,1),...
            'muTimeStamp',NaN(NBblock,1),...
            'altimeter',NaN(NBblock,2),...
            'fsync_err',NaN(NBblock,1),...
            'Checksum_aux1',NaN(NBblock,1),...
            'Checksum_map',NaN(NBblock,1));
    otherwise
        EPSI.madre = struct(...
            'EpsiStamp',NaN(NBblock,1),...
            'TimeStamp',NaN(NBblock,1),...
            'altimeter',NaN(NBblock,2),...
            'fsync_err',NaN(NBblock,1),...
            'Checksum_aux1',NaN(NBblock,1),...
            'Checksum_map',NaN(NBblock,1));
end

if(isfield(EPSI,'header'))
    EPSI.madre.time = NaN(NBblock,1);
end

if is_aux1
    EPSI.aux1 = struct(...
        'Aux1Stamp',NaN(NBblock*aux1.nbsample,1),...
        'T_raw',NaN(NBblock*aux1.nbsample,1),...
        'C_raw',NaN(NBblock*aux1.nbsample,1),...
        'P_raw',NaN(NBblock*aux1.nbsample,1),...
        'PT_raw',NaN(NBblock*aux1.nbsample,1));
end

% for cha=1:Meta_Data.PROCESS.nb_channels
%     wh_channel=Meta_Data.PROCESS.channels{cha};
%     EPSI.epsi.(wh_channel) = NaN(NBblock,epsi.nbsamples);
% end
EPSI.epsi.EPSInbsample=NaN(NBblock,epsi.nbsamples);

% we are checking if the very last block is good too.
check_endstr=mod(ind_madre(end)+epsi.offset(end)+epsi.name_length+epsi.total_length ...
                                - numel(str),epsi.total_length);
if check_endstr==0
    nb_block=NBblock;
else
    nb_block=NBblock-1;
end
% still initializing
if is_aux1
    aux1.stamp = char(zeros(nb_block*aux1.nbsample,aux1.stamp_length));
    aux1.sbe = char(zeros(nb_block*aux1.nbsample,aux1.sbe_length));
end


datalength=epsi.nbsamples*epsi.nchannels*epsi.bytes_per_adc; % length of an EPSI block with 3 bytes. For some reason the epsi data are in hex and mess up the byte count
epsi.raw = int32(zeros(nb_block,datalength));


% now lets begin  reading and splitting!!!!
for i=1:numel(indblock)
    % grab local time if STREAMING SITUATION;
    if(isfield(EPSI,'header'))
        if (isfield(EPSI.header,'offset_time'))
            system.time(i,1:10) = str(ind_madre(i)-(10:-1:1));
        end
    end
    % read items in the EPSI block Header
    madre.epsi_stamp(i,:) = str(indblock(i)+(1:madre.epsi_stamp_length)+madre.epsi_stamp_offset);
    madre.epsi_time(i,:) = str(indblock(i)+(1:madre.epsi_time_length)+madre.epsi_time_offset);
    madre.aux1_chksum(i,:) = str(indblock(i)+(1:madre.aux_chksum_length)+madre.aux_chksum_offset);
    madre.epsi_chksum(i,:) = str(indblock(i)+(1:madre.map_chksum_length)+madre.map_chksum_offset);
    for j=1:2
        madre.altimeter((i-1)*2+j,:) = str(indblock(i)+(1:madre.alt_time_length)+madre.alt_time_offset(j));
    end
    madre.fsync_err(i,:) = str(indblock(i)+(1:madre.fsync_err_length)+madre.fsync_err_offset);
    switch firmware_version
        case 'microsecond' % ALB: I added the musecond timestamp in third position. So the header is longer
            madre.epsi_mutime(i,:) = str(indblock(i)+(1:madre.epsi_time_length)+madre.epsi_mutime_offset);
    end
    % read AUX1 (SBE49) block 
    if is_aux1
        for j=1:aux1.nbsample
            aux1.stamp((i-1)*aux1.nbsample+j,:) = str(indblock(i)+(1:aux1.stamp_length)+aux1.stamp_offset(j));
            aux1.sbe((i-1)*aux1.nbsample+j,:) = str(indblock(i)+(1:aux1.sbe_length)+aux1.sbe_offset(j));
        end
    end
    % get the EPSI block
    if hexepsi==0
        fprintf("block #%i over %i blocks \r\n",i,numel(indblock))
        epsi.raw(i,:) = int32(str(indblock(i)+epsi.offset(end)+epsi.name_length+(1:epsi.total_length)));
        try
            EPSI.aux1.T_raw(i) = hex2dec(aux1.sbe(i,1:6));
        catch
            disp('issue with hex')
            idxstamp=aux1.sbe(i,1:6)<0; % make a 0 arrays
            hexidx=regexp(aux1.sbe(i,1:6),'([0-9A-Fa-f])');
            idxstamp(hexidx)=1;
            aux1.sbe(i,idxstamp==0)= ...
                aux1.sbe(i-1,idxstamp==0);
            EPSI.aux1.T_raw(i) = hex2dec(aux1.sbe(i,1:6));
            
        end
    else
        fprintf("block #%i over %i blocks \r\n",i,numel(indblock))
        tempo_str=str(indblock(i)+epsi.offset(end)+epsi.name_length+(1:epsi.total_length));
        tempo1=strsplit(tempo_str,'\r\n');
        test=[tempo1{1:end-1}];
        count=0;
        for j=1:2:length(test)
            count=count+1;
            epsi.raw(i,count) = hex2dec(test(j:j+1));
        end
    end

end
% done with split file
%toc

%convert 3 bytes ADC samples into 24 bits counts. 



epsi.raw1 = epsi.raw(:,1:epsi.bytes_per_adc:end)*256^2+ ...
            epsi.raw(:,2:epsi.bytes_per_adc:end)*256+ ...
            epsi.raw(:,3:epsi.bytes_per_adc:end);


if(isfield(EPSI,'header'))
    switch Meta_Data.PROCESS.recording_mode
        case 'STREAMING'
            system.time(:,11) = newline;
            system.time = system.time';
            system_time = textscan(system.time(:),'%f');
            EPSI.madre.time = system_time{1}/100/24/3600+EPSI.header.offset_time;
        case 'SD'
            EPSI.madre.time=0;
    end
end

% converting Hex into decimal. Starts with the header.
try
    EPSI.madre.EpsiStamp = hex2dec(madre.epsi_stamp);
    EPSI.madre.TimeStamp = hex2dec(madre.epsi_time);
    EPSI.madre.altimeter = reshape(hex2dec(madre.altimeter),2,[])';
    EPSI.madre.fsync_err = hex2dec(madre.fsync_err);
    EPSI.madre.Checksum_aux1 = hex2dec(madre.aux1_chksum);
    EPSI.madre.Checksum_map = hex2dec(madre.epsi_chksum);
    switch firmware_version
        case 'microsecond' % ALB: I added the musecond timestamp in third position. So the header is longer
            EPSI.madre.muTimeStamp = hex2dec(madre.epsi_mutime);
    end
catch
    disp('issue with madre time')
    mask=madre.epsi_stamp==' ';madre.epsi_stamp(mask)='0';
    for ii=1:length(indblock)
        try
            EPSI.madre.EpsiStamp(ii)=...
                hex2dec(madre.epsi_stamp(ii,:));
        catch
            disp('toto')
            idxstamp=madre.epsi_stamp(ii,:)<0; % make a 0 arrays 
            hexidx=regexp(madre.epsi_stamp(ii,:),'([0-9A-Fa-f])');
            idxstamp(hexidx)=1;
            madre.epsi_stamp(ii,idxstamp==0)=...
                madre.epsi_stamp(ii-1,idxstamp==0);
            EPSI.madre.EpsiStamp(ii)=...
                hex2dec(madre.epsi_stamp(ii,:));
        end
        try
            EPSI.madre.TimeStamp(ii)=...
                hex2dec(madre.epsi_time(ii,:));
        catch
            disp('tata')
            idxstamp=madre.epsi_time(ii,:)<0; % make a 0 arrays 
            hexidx=regexp(madre.epsi_time(ii,:),'([0-9A-Fa-f])');
            idxstamp(hexidx)=1;
            madre.epsi_time(ii,idxstamp==0)= ...
                madre.epsi_time(ii-1,idxstamp==0);
            EPSI.madre.TimeStamp(ii)=...
                hex2dec(madre.epsi_time(ii,:));
        end

    end
%     EPSI.madre.EpsiStamp = hex2dec(madre.epsi_stamp);
%     EPSI.madre.TimeStamp = hex2dec(madre.epsi_time);
%     EPSI.madre.altimeter = reshape(hex2dec(madre.altimeter),2,[])';
%     EPSI.madre.fsync_err = hex2dec(madre.fsync_err);
%     EPSI.madre.Checksum_aux1 = hex2dec(madre.aux1_chksum);
%     EPSI.madre.Checksum_map = hex2dec(madre.epsi_chksum);
end
% issues with the SD write and some bytes are not hex. if issues we scan
% the whole sbe time series to find the bad bytes and then use the average 
% increment from with the previous samples;
% TO DO get rid of the nameam Tdiff over 10s.
if is_aux1
    try 
        EPSI.aux1.T_raw = hex2dec(aux1.sbe(:,1:6));
        EPSI.aux1.C_raw = hex2dec(aux1.sbe(:,(1:6)+6));
        EPSI.aux1.P_raw = hex2dec(aux1.sbe(:,(1:6)+12));
        EPSI.aux1.PT_raw = hex2dec(aux1.sbe(:,(1:4)+18));
        EPSI.aux1.Aux1Stamp =hex2dec(aux1.stamp);
    catch 
        disp('bug in SBE hex bytes')
        % ALB:San's trick to get the errors and replace the bad char by '0'
        ind_sbe = ( aux1.sbe >='0' & ...
                aux1.sbe <='9')| ...
              ( aux1.sbe >='a' & ...
                aux1.sbe <='f')| ...
              ( aux1.sbe >='A' & ...
                aux1.sbe <='F');
            
        ind_stamp = ( aux1.stamp >='0' & ...
                aux1.stamp <='9')| ...
              ( aux1.stamp >='a' & ...
                aux1.stamp <='f')| ...
              ( aux1.stamp >='A' & ...
                aux1.stamp <='F');
            % replace bad char with '0'
            aux1.stamp(~ind_stamp)='0';
            aux1.sbe(~ind_sbe)='0';
            EPSI.aux1.T_raw = hex2dec(aux1.sbe(:,1:6));
            EPSI.aux1.C_raw = hex2dec(aux1.sbe(:,(1:6)+6));
            EPSI.aux1.P_raw = hex2dec(aux1.sbe(:,(1:6)+12));
            EPSI.aux1.PT_raw = hex2dec(aux1.sbe(:,(1:4)+18));
            EPSI.aux1.Aux1Stamp =hex2dec(aux1.stamp);
            
    end

    [EPSI.aux1.Aux1Stamp,ia0,~] =unique(EPSI.aux1.Aux1Stamp,'stable');
    EPSI.aux1.Aux1Stamp=filloutliers(EPSI.aux1.Aux1Stamp,'center','movmedian',100);
    %ALB reorder the stamps and samples because until now we kept the zeros
    % in the aux block
    [EPSI.aux1.Aux1Stamp,ia1]=sort(EPSI.aux1.Aux1Stamp);
    EPSI.aux1.T_raw  = EPSI.aux1.T_raw(ia0(ia1));
    EPSI.aux1.C_raw  = EPSI.aux1.C_raw(ia0(ia1));
    EPSI.aux1.P_raw  = EPSI.aux1.P_raw(ia0(ia1));
    EPSI.aux1.PT_raw = EPSI.aux1.PT_raw(ia0(ia1));
    
    EPSI = epsi_ascii_get_temperature(EPSI);
    EPSI = epsi_ascii_get_pressure(EPSI);
    EPSI = epsi_ascii_get_conductivity(EPSI);

    % remove bad records for aux1
    ind = EPSI.aux1.Aux1Stamp == 0 & EPSI.aux1.T_raw == 0 & EPSI.aux1.C_raw == 0 & EPSI.aux1.P_raw == 0;
    aux1_fields = fieldnames(EPSI.aux1);
    for i  = 1:numel(aux1_fields)
        EPSI.aux1.(aux1_fields{i})(ind) = NaN;
    end
end


% parsing the EPSI block data
for cha=1:Meta_Data.PROCESS.nb_channels
    wh_channel=Meta_Data.PROCESS.channels{cha};
    EPSI.epsi.([wh_channel '_count']) = epsi.raw1(:,cha:epsi.nchannels:end);
end

% input epsi sample stamp on each according to the sample sent via madre
% record
EPSI.epsi.EPSInbsample = repmat(1:epsi.nbsamples,[NBblock 1])+repmat(EPSI.madre.EpsiStamp,[1 epsi.nbsamples])-epsi.nbsamples;
if(isfield(EPSI,'header'))
    % delayed by 1 sample period assuming that it would take a bit of time
    % to transfer the data.
    switch Meta_Data.PROCESS.recording_mode
        case 'STREAMING'
            EPSI.epsi.time = (repmat(1:epsi.nbsamples,[NBblock 1])-epsi.nbsamples-1)*epsi.sample_period/24/3600 + repmat(EPSI.madre.time,[1 epsi.nbsamples]);
        case 'SD'
            EPSI.epsi.time = Meta_Data.starttime+EPSI.epsi.EPSInbsample/Meta_Data.PROCESS.Fs_epsi/86400;
    end
end

% coefficient for Unipolar or bipolar ADC configuration and also to convert
% accelerometer Voltage into Accelereation units (in g).
full_range = 2.5;
bit_counts = 24;
gain = 1;
acc_offset = 1.65;
acc_factor = 0.66;

for cha=1:Meta_Data.PROCESS.nb_channels
    wh_channel=Meta_Data.PROCESS.channels{cha};
    %if ~strcmp(wh_channel,'c') %NC 10/12/21 - change to 'c_count')
    if ~strcmp(wh_channel,'c_count')
        switch Meta_Data.epsi.(wh_channel).ADCconf
            case {'Bipolar','bipolar'}
                EPSI.epsi.([wh_channel '_volt'])=full_range/gain* ...
                    (double(EPSI.epsi.([wh_channel '_count']))/2.^(bit_counts-1)-1);
            case {'Unipolar','unipolar'}
                EPSI.epsi.([wh_channel '_volt'])=full_range/gain* ...
                    double(EPSI.epsi.([wh_channel '_count']))/2.^(bit_counts);
                
        end
        
        switch wh_channel
            case 'a1'
                EPSI.epsi.a1_g = (EPSI.epsi.a1_volt-acc_offset)/acc_factor;
            case 'a2'
                EPSI.epsi.a2_g = (EPSI.epsi.a2_volt-acc_offset)/acc_factor;
            case 'a3'
                EPSI.epsi.a3_g = (EPSI.epsi.a3_volt-acc_offset)/acc_factor;
        end
    end
end
% Remove the acceleration volts fields
try
    EPSI.epsi = rmfield(EPSI.epsi,'a1_volt');
catch
    disp("no a1_volt")
end
try
EPSI.epsi = rmfield(EPSI.epsi,'a2_volt');
catch
    disp("no a2_volt")
end
try
EPSI.epsi = rmfield(EPSI.epsi,'a3_volt');
catch
    disp("no a3_volt")
end

% grab all the epsi field names.
epsi_fields = fieldnames(EPSI.epsi);


% lay all records out straight instead of bunching them up with the MADRE
% records
for i  = 1:numel(epsi_fields)
    EPSI.epsi.(epsi_fields{i}) = reshape(EPSI.epsi.(epsi_fields{i})',[],1);
end

if(isfield(EPSI.epsi,'time') && ~isempty(EPSI.epsi.time)) && is_aux1
    EPSI.aux1.time = NaN(size(EPSI.aux1.Aux1Stamp));
    [ia, ib] = ismember(EPSI.aux1.Aux1Stamp,EPSI.epsi.EPSInbsample);
    EPSI.aux1.time(ia) = EPSI.epsi.time(ib(ib>0));
end
%toc
end


%%
function epsi = mod_combine_epsi(varargin)
% mod_combine_epsi - combines epsi data files in MATLAB format that was
% converted using mod_read_epsi_raw
%
% mod_combine_epsi(epsi1,epsi2,epsi3,...) returns a EPSI structure of variables described
% for MET data files
%
% Written 2018/10/15 - San Nguyen stn 004@ucsd.edu


if nargin < 1
    epsi = [];
    return;
end

if nargin == 1
    epsi = varargin{1};
    if length(epsi) == 1
        return;
    end
    evalstr = 'epsi = mod_combine_epsi(';
    for i=1:(length(epsi)-1)
        evalstr = strcat(evalstr, 'epsi(', num2str(i), '),');
    end
    evalstr = strcat(evalstr, 'epsi(', num2str(length(epsi)), '));');
    eval(evalstr);
    return
end

epsi_fields = fieldnames(varargin{1});
for i = 2:nargin
    tmp_fields = fieldnames(varargin{i});
    for j = 1:numel(tmp_fields)
        if ~ismember(tmp_fields{j},epsi_fields)
            epsi_fields{end+1} = tmp_fields{j};
        end
    end
end

epsi_sub_fields = cell(size(epsi_fields));

for i = 1:numel(epsi_fields)
    epsi_sub_fields{i} = fieldnames(varargin{1}.(epsi_fields{i}));
    
    for j = 2:nargin
        tmp_fields = fieldnames(varargin{j}.(epsi_fields{i}));
        for k = 1:numel(tmp_fields)
            if ~ismember(tmp_fields{k},epsi_sub_fields{i})
                epsi_sub_fields{i}{end+1} = tmp_fields{k};
            end
        end
    end
end

for i=1:(length(epsi_fields))
    %header field
    if strcmpi(epsi_fields{i},'header')
        evalstr = strcat('epsi.', epsi_fields{i}, '= [');
        for j=1:(nargin-1)
            if ~isfield(varargin{j},(epsi_fields{i}))
                varargin{j}.(epsi_fields{i}) = NaN(size(varargin{j}.Time));
            end
            evalstr = strcat(evalstr, 'varargin{', num2str(j), '}.', epsi_fields{i}, ';');
        end
        if ~isfield(varargin{nargin},(epsi_fields{i}))
            varargin{nargin}.(epsi_fields{i}) = NaN(size(varargin{nargin}.Time));
        end
        evalstr = strcat(evalstr, 'varargin{', num2str(nargin), '}.', epsi_fields{i}, '];');
        eval(evalstr);
        continue;
    end
    % other fields
    for j=1:(length(epsi_sub_fields{i}))
        evalstr = strcat('epsi.',epsi_fields{i} ,'.', epsi_sub_fields{i}{j}, '= [');
        for k=1:(nargin-1)
            if ~isfield(varargin{k},epsi_fields{i})
                continue;
            end
            if ~isfield(varargin{k}.(epsi_fields{i}),(epsi_sub_fields{i}{j}))
                continue;
            end
            evalstr = strcat(evalstr, 'varargin{', num2str(k), '}.',epsi_fields{i} ,'.',epsi_sub_fields{i}{j}, ';');
        end
        if ~isfield(varargin{nargin},(epsi_fields{i}))
            evalstr = strcat(evalstr, '];');
        elseif ~isfield(varargin{nargin}.(epsi_fields{i}),(epsi_sub_fields{i}{j}))
            evalstr = strcat(evalstr, '];');
        else
            evalstr = strcat(evalstr, 'varargin{', num2str(nargin), '}.',epsi_fields{i} ,'.', epsi_sub_fields{i}{j}, '];');
        end
        eval(evalstr);
    end
end
end

%  parse all the lines in the header of the file
function EPSI = epsi_ascii_parseheader(fid)
EPSI = [];
fgetl(fid);
s=fgetl(fid);
[v,val]=epsi_ascii_parseheadline(s);
if ~isempty(v)
    eval(['EPSI.header.' lower(v) '=' val ';']);
end
s=fgetl(fid);
if(~strncmp(s,'%*****START_FCTD',16))
    return;
end

s=fgetl(fid);
while ~strncmp(s,'%*****END_FCTD',14) && ~feof(fid)
    [v,val]=epsi_ascii_parseheadline(s);
    if ~isempty(v)
        try
            eval(['EPSI.header.' lower(v) '=' val ';']);
        catch obj
            if strncmp(v,'FCTD_VER',8)
                eval(['EPSI.header.' lower(v) '=''' val ''';']);
            else
                %                 disp(obj.message);
                %                 disp(['Error occured in string: ' s]);
            end
            
        end
    end
    s=fgetl(fid);
    %     strncmp(s,'%*****END_FCTD',14);
end
return;
end

%  parse each line in the header to detect comments
function [v,val]=epsi_ascii_parseheadline(s)
if(isempty(s))
    v = [];
    val = [];
    return;
end
if s(1)~='%'
    
    i = strfind(s,'=');
    v=s(1:i-1);
    val = s(i+1:end);
else
    v=[];
    val=[];
end

return;
end


function corrected_time = epsi_ascii_correct_negative_time(time)
corrected_time = time;
neg_time = time(time<0);
corrected_time(time<=0) = 2^64 + neg_time;
end


%  reads and apply calibration to the temperature data
function EPSI = epsi_ascii_get_temperature(EPSI)

a0 = EPSI.header.ta0;
a1 = EPSI.header.ta1;
a2 = EPSI.header.ta2;
a3 = EPSI.header.ta3;

mv = (EPSI.aux1.T_raw-524288)/1.6e7;
r = (mv*2.295e10 + 9.216e8)./(6.144e4-mv*5.3e5);
EPSI.aux1.T = a0+a1*log(r)+a2*log(r).^2+a3*log(r).^3;
EPSI.aux1.T = 1./EPSI.aux1.T - 273.15;
EPSI.aux1.T=real(EPSI.aux1.T);
return;
end

%  reads and apply calibration to the conductivity data
function EPSI = epsi_ascii_get_conductivity(EPSI)
try 
g = EPSI.header.g;
h = EPSI.header.h;
i = EPSI.header.i;
j = EPSI.header.j;
tcor = EPSI.header.tcor;
pcor = EPSI.header.pcor;
catch
g = EPSI.header.cg;
h = EPSI.header.ch;
i = EPSI.header.ci;
j = EPSI.header.cj;
tcor = EPSI.header.ctcor;
pcor = EPSI.header.cpcor;
end

f = EPSI.aux1.C_raw/256/1000;

EPSI.aux1.C = (g+h*f.^2+i*f.^3+j*f.^4)./(1+tcor*EPSI.aux1.T+pcor*EPSI.aux1.P);
EPSI.aux1.C=real(EPSI.aux1.C);
return;
end

%  reads and apply calibration to the pressure data
function EPSI = epsi_ascii_get_pressure(EPSI)
% ALB 04112019 Changed EPSI.header.SBEcal. to EPSI.header.
pa0 = EPSI.header.pa0;
pa1 = EPSI.header.pa1;
pa2 = EPSI.header.pa2;
ptempa0 = EPSI.header.ptempa0;
ptempa1 = EPSI.header.ptempa1;
ptempa2 = EPSI.header.ptempa2;
ptca0 = EPSI.header.ptca0;
ptca1 = EPSI.header.ptca1;
ptca2 = EPSI.header.ptca2;
ptcb0 = EPSI.header.ptcb0;
ptcb1 = EPSI.header.ptcb1;
ptcb2 = EPSI.header.ptcb2;


y = EPSI.aux1.PT_raw/13107;

t = ptempa0+ptempa1*y+ptempa2*y.^2;
x = EPSI.aux1.P_raw-ptca0-ptca1*t-ptca2*t.^2;
n = x*ptcb0./(ptcb0+ptcb1*t+ptcb2*t.^2);

EPSI.aux1.P = (pa0+pa1*n+pa2*n.^2-14.7)*0.689476;

return;
end