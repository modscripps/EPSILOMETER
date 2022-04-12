function [epsi,ctd,alt,act,vnav,gps,seg,spec] = mod_som_read_epsi_files_v3(filename,Meta_Data)

% For now, only works for 'v3' of som_acq
% It only works with a single file, since we're always calling it from a
% loop  in Epsi_MakeMatFromRaw
%
% TO DO:
%   Add backwards compatability - look at mod_som_read_epsi_files_v2.m
%
% CALLED BY:
%   Epsi_MakeMatFromRaw
%
% INPUTS
%   filename
%   Meta_Data
%
% OUTPUTS
%   epsi  (intermediate processing in 'efe' structure)
%   ctd   (intermediate processing in 'sbe' structure)
%   alt   (intermediate processing in 'alti' structure)
%   act   (intermediate processing in 'actu' structure)
%   vnav  (intermediate processing in 'vecnav' structure)
%   gps   (intermediate processing in 'gpsmeta' structure
%
% Nicole Couto adapted from Arnaud LeBoyer's mod_som_read_epsi_raw.m
% June 2021
% -------------------------------------------------------------------------

% Define a constant for salinity calculation
c3515 = 42.914;

%% Open file and save contents as 'str'
fid = fopen(filename);
fseek(fid,0,1);
frewind(fid);
str = fread(fid,'*char')';
fclose(fid);

% To do: Add the other versions and a switch to choose between them later. For
% now, use 'v3' only.
efe_block_version = 'v3';

%% Get indices and tokens for each data type you will process
% ind_*_start  = starting indices of all matches
% ind_*_end    = ending indices of all matches

[ind_efe_start, ind_efe_stop]   = regexp(str,'\$EFE([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
[ind_sbe_start, ind_sbe_stop]   = regexp(str,'\$SB49([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
[ind_alt_start,ind_alt_stop]    = regexp(str,'\$ALT([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
[ind_act_start,ind_act_stop]    = regexp(str,'\$ACTU([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
[ind_vnav_start,ind_vnav_stop]  = regexp(str,'\$VNMAR([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
[ind_gps_start, ind_gps_stop]   = regexp(str,'\$GPGGA([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');

[ind_seg_start, ind_seg_stop]   = regexp(str,'\$SEGM([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
[ind_spec_start, ind_spec_stop]         = regexp(str,'\$SPEC([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
[ind_avgspec_start, ind_avgspec_stop]   = regexp(str,'\$AVGS([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');
[ind_dissrate_start, ind_dissrate_stop] = regexp(str,'\$RATE([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end');


    %ALB GROSCON
    spltfilename=strsplit(filename,'/');
    strdate=spltfilename{end}(12:end-4);
    strdate=['20' strdate];
    t0=(datenum(strdate,"yyyy_mm_dd_HHMMSS")-datenum("01-01-1970 00:00:00"))*86400;
    t0=0;

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



%% Process EFE data
if isempty(ind_efe_start)
    disp('no epsi data')
    epsi=[];
else
    disp('processing epsi data')
    % EFE-specific quantities
    efe.data.n_blocks           = numel(ind_efe_start); %Number of data blocks beginning with $EFE
    
    efe.data.n_channels         = 7;
    efe.data.sample_freq        = 320;
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
    % epsi.time_s and epsi.dnum will be created from epsi_timestamp once
    % all its records are filled
    epsi_timestamp   = nan(efe.data.n_recs,1);
   
    
    
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
    %ALB GROSCON
%     epsi_timestamp=epsi_timestamp-epsi_timestamp(1)+t0*1000;
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
    epsi = orderfields(epsi,{'dnum','time_s','t1_count','t2_count','s1_count',...
        's2_count','a1_count','a2_count','a3_count','t1_volt','t2_volt','s1_volt',...
        's2_volt','a1_g','a2_g','a3_g'});
    
end %end if there is epsi data


%% CTD
if isempty(ind_sbe_start)
    disp('no ctd data')
    ctd=[];
else
    disp('processing ctd data')
    
    % SBE-specific quantities
    % ---------------------------------
    
    switch Meta_Data.CTD.name
        case{"SBE49","SBE","S49","SB49"}
            sbe.data.format      = 'eng';
            sbe.data.length      = 22;
            sbe.data.sample_freq = 16;
            sbe.cal              = Meta_Data.CTD.cal;
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
        if (length(sbe_block_data)~=(24+16)*sbe.data.recs_per_block)
            fprintf("SBE block %i has incorrect length\r\n",iB)
        else
            %
            % Where do 16 and 24 come from?
            sbe_block_data=reshape(sbe_block_data,16+24,sbe.data.recs_per_block).';
            
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
                        ctd.P(n_rec)      = str2double(data_split{1});
                        ctd.T(n_rec)      = str2double(data_split{2});
                        ctd.S(n_rec)      = str2double(data_split{3});
                        ctd.C(n_rec)      = NaN;
                    case 'eng'
                        ctd.T_raw(n_rec)  = hex2dec(rec_ctd(:,1:6));
                        ctd.C_raw(n_rec)  = hex2dec(rec_ctd(:,(1:6)+6));
                        ctd.P_raw(n_rec)  = hex2dec(rec_ctd(:,(1:6)+12));
                        ctd.PT_raw(n_rec) = hex2dec(rec_ctd(:,(1:4)+18));
                end %end switch sbe.data.format
            end %end loop through recs_per_block
        end %end if sbe data block is the correct size
    end %end loop through sbe blocks
    
    % If timestamp has values like 1.6e12, it is in milliseconds since Jan
    % 1, 1970. Otherwise it's in milliseconds since the start of the record
        %ALB GROSCON
%     ctd_timestamp=ctd_timestamp-ctd_timestamp(1)+t0*1000;

    if nanmedian(ctd_timestamp)>1e9
        % time_s - seconds since 1970
        % dnum - matlab datenum
        [ctd.time_s,ctd.dnum] = convert_timestamp(ctd_timestamp);
        disp(datestr(ctd.dnum(1)))
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
            ctd.th   = sw_ptmp(ctd.S,ctd.T,ctd.P,0);
            ctd.sgth  = sw_pden(ctd.S,ctd.T,ctd.P,0);
            ctd.dPdt = [0; diff(ctd.P)./diff(ctd.time_s)];
            
            % NC 17 July 2021 - added ctd.z and ctd.dzdt.
            % get_scan_spectra.m will use dzdt to define fall speed w.
            if ~isfield(Meta_Data.PROCESS,'latitude')
                error('Need latitude to get depth from pressure data. Add to MetaProcess text file.')
            else
                ctd.z    = sw_dpth(ctd.P,Meta_Data.PROCESS.latitude);
                ctd.dzdt = [0; diff(ctd.z)./diff(ctd.time_s)];
            end

    end
    
    % Make sure ctd.S is real. Every once in a while, S comes out imaginary, I think only when SBE is on deck.
    ctd.S = real(ctd.S);
    
    % Sort ctd fields
    ctd = orderfields(ctd,{'dnum','time_s','P_raw','T_raw','S_raw',...
        'C_raw','PT_raw','P','z','T','S','C','th','sgth','dPdt','dzdt'});

end %end loop if there is ctd data

%% Altimeter
if isempty(ind_alt_start)
    disp('no alt data')
    alt=[];
else
    disp('processing alt data')
    
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
        alti.hextimestamp.value   = hex2dec(alt_block_str(tag.hextimestamp.inds));
        alt_timestamp(iB) = alti.hextimestamp.value;

        
        % The altimeter block does not have a hexlengthblock (hexadecimal
        % length of data block). It has the altimeter reading in that
        % position
        alti.hexlengthblock.value = alt_block_str(tag.hexlengthblock.inds);
        dst_microsec = str2double(alti.hexlengthblock.value);
        alt.dst(iB) = dst_microsec*alti.ten_microsec2sec*alti.sound_speed;
        
    end %end loop through alt blocks
    
    % If timestamp has values like 1.6e12, it is in milliseconds since Jan
    % 1, 1970. Otherwise it's in milliseconds since the start of the record
    %GROSCON
%     alt_timestamp=alt_timestamp-alt_timestamp(1)+t0*1000;
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
    alt = orderfields(alt,{'dnum','time_s','dst'});
    
end %end loop if there is alt data


%% Actuator
if isempty(ind_act_start)
    disp('no act data')
    act=[];
else
    disp('processing act data')
    
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
    disp('no vnav data')
    vnav=[];
else
    disp('processing vnav data')
    
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
    vnav.compass = nan(vecnav.data.n_recs,1);
    vnav.acceleration = nan(vecnav.data.n_recs,1);
    vnav.gyro = nan(vecnav.data.n_recs,1);
    
%     % Grab the block of data starting with the header
%     vnav_block_str = str(ind_vnav_start(iB):ind_vnav_stop(iB));
%     %Commented out and moved below by Bethan June 26
    
    
    % Loop through data blocks and parse strings
    for iB=1:vecnav.data.n_blocks
        
        % Grab the block of data starting with the header
        vnav_block_str = str(ind_vnav_start(iB):ind_vnav_stop(iB)); %Moved here by Bethan June 26
        % For the vecnav, we use the hexadecimal timestamp that comes before
        % $VNMAR
        vnav_timestamp(iB) = hex2dec(str(ind_time_start(iB):ind_time_stop(iB)));

        % Get the data after the header
        vnav_block_data = str(ind_vnav_start(iB):ind_vnav_stop(iB)-tag.chksum.length);
        
        % Split the data into parts
        data_split = strsplit(vnav_block_data,',');
        
        % Compass, acceleration, and gyro each have x, y, and z components
        for iC=1:3
            vnav.compass(iB,iC)=str2double(data_split{1,iC+1}); %Compass (x,y,z) units = gauss
            vnav.acceleration(iB,iC)=str2double(data_split{1,iC+4}); %Acceleration (x,y,x) units = m/s^2
            vnav.gyro(iB,iC)=str2double(data_split{1,iC+7}); %Gyro (x,y,z) units = rad/s
        end
        
    end
    
    % If timestamp has values like 1.6e12, it is in milliseconds since Jan
    % 1, 1970. Otherwise it's in milliseconds since the start of the record
%     vnav_timestamp=vnav_timestamp-vnav_timestamp(1)+t0*1000;
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
    vnav = orderfields(vnav,{'dnum','time_s','compass','acceleration','gyro'});
    
end %end loop if there is vnav data

%% GPS data
if isempty(ind_gps_start)
    disp('no gps data')
    gps = [];
else
    
    disp('processing gps data')
    
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
    
    try
    % Loop through data blocks and parse strings
    for iB=1:gpsmeta.data.n_blocks
        
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
        
    end
    
    catch
            warning("gps issue")
        gps.dnum(iB)      = nan;
        gps.latitude(iB)  = nan;
        gps.longitude(iB) = nan;
        gps.longitude(iB) = nan;
        end
    
end %end loop if there is gps data

%% Process SEGM data
if isempty(ind_seg_start)
    disp('no segment data')
    seg=[];
else
    disp('processing segment data')
    
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
        seg.time_s = arrayfun(@(x,y) linspace(x,y, ...
                                             seg.data.sample_per_segment)./1000, ...
                                    segment_timestamp(1:end-1), ...
                                    segment_timestamp(2:end),'un',0);
                                
        % ALB 500 is 1000*0.5. 
        % ALB 1000 because segment_timestampare in milisec/
        % ALB 0.5 becasue I have 50% overlap.
        seg.time_s{end+1}=linspace(segment_timestamp(end), ...
                                   segment_timestamp(end)+500*seg.data.sample_per_segment/seg.data.sample_freq,...
                                   seg.data.sample_per_segment)./1000;
%         seg.time_s=cell2mat(seg.time_s.');
         seg.dnum=cellfun(@(x) x.*nan,seg.time_s,'un',0);
    end
    
    % Sort epsi fields
    seg = orderfields(seg,{'dnum','time_s','t1_volt','s1_volt','a3_g','data','channels','checksum'});
    
end %end if there is epsi data

%% Process SPEC data
if isempty(ind_spec_start)
    disp('no spec data')
    spec=[];
else
    disp('processing spec data')
    
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
    
    if nanmedian(spec_timestamp)>1e9
        % time_s - seconds since 1970
        % dnum - matlab datenum
        [spec.time_s,spec.dnum] = convert_timestamp(spec_timestamp);
    else
        % time_s - seconds since power on
        spec.time_s = spec_timestamp./1000;
        spec.dnum=spec.time_s.*nan;

    end
    
    % Sort epsi fields
    spec = orderfields(spec,{'dnum','time_s','t1_volt','s1_volt','a3_g','data','channels','checksum'});
    
end %end spec






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
                disp(obj.message);
                disp(['Error occured in string: ' s]);
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

