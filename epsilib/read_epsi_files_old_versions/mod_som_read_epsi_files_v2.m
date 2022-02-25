function [epsi,ctd,alt,act]=mod_som_read_epsi_files_v2(filename,Meta_Data)

% check if it is a single file or a directory and a set of files
if ischar(filename) % dir or file
    switch exist(filename,'file')
        case 2 % if it is a file
            %             fid = fopen(filename,'r');
            [epsi,ctd,alt,act] = mod_som_read_epsi_raw(filename,Meta_Data);
            %             fclose(fid);
        case 7 % if it is a directory
            fileDir = filename;
            
            %fileList = [];
            %fileList = [fileList; dir(fullfile(fileDir,'*_raw'))];
            fileList = [fileList; dir(fullfile(fileDir,'*.raw'))];
            %fileList = [fileList; dir(fullfile(fileDir,['*',Meta_Data.rawfileSuffix]))];

            if isempty(fileList)
                epsi = [];
                return
            else
                % prepare to read all files
                epsi = cell(size(fileList));
                ctd  = cell(size(fileList));
                alt  = cell(size(fileList));
                act = cell(size(fileList));
                vnav = cell(size(fileList));
                % read the files in the directory
                for i = 1:length(fileList)
                    
                    % Read data in current file
                    disp(['reading ' fileList(i).name]);
                    [epsi{i},ctd{i},alt{i},act{i}] = mod_som_read_epsi_raw(fullfile(fileDir,fileList(i).name),Meta_Data);
                    
                    structList = {'epsi','ctd','alt','vnav'};
                    %NC - add tracker for file
                    %for iS=1:length(structList)
                    %    if isstruct(eval(structList{iS}))
                    %end
                    
                    if isstruct(epsi{i})
                        epsi{i}.fileNum = i*ones(size(epsi{i}.timestamp,1),size(epsi{i}.timestamp,2));
                    end
                    if isstruct(ctd{i})
                        ctd{i}.fileNum = i*ones(size(ctd{i}.timestamp,1),size(ctd{i}.timestamp,2));
                    end
                    if isstruct(alt{i})
                        alt{i}.fileNum = i*ones(size(alt{i}.timestamp,1),size(alt{i}.timestamp,2));
                    end
                    
                end
                % combine all files into one structure
                if sum(cellfun(@length,epsi))>0
                    epsi = mod_combine_epsi(epsi{:});
                end
                if sum(cellfun(@length,ctd))>0
                    ctd = mod_combine_ctd(ctd{:});
                end
                if sum(cellfun(@length,alt))>0
                    alt = mod_combine_alt(alt{:});
                end
            end
        otherwise
            error('MATLAB:mod_read_epsi_raw:wrongFILENAME','Invalid file specified.');
    end
elseif iscellstr(filename) % cell of files
    % prepare to read all files
    epsi = cell(size(filename));
    ctd  = cell(size(filename));
    alt  = cell(size(filename));
    % read all files
    for i = 1:length(filename)
        disp(['reading ' filename{i}]);
        [epsi{i},ctd{i},alt{i}] = mod_som_read_epsi_raw(filename{i},Meta_Data);
    end
    % combine all files into one epsi structure
    epsi = mod_combine_epsi(epsi{:});
    ctd  = mod_combine_ctd(ctd{:});
    alt  = mod_combine_alt(alt{:});
else
    
    if (filename<1)
        error('MATLAB:mod_read_epsi_raw:wrongFID','FID is invalid');
    end
    
    [epsi,ctd,alt] = mod_som_read_epsi_raw(filename,Meta_Data);
    
    return
end

% save(fullfile(Meta_Data.Epsipath,['epsi_' Meta_Data.deployment '.mat']),'epsi','-v7.3');
% save(fullfile(Meta_Data.CTDpath, ['ctd_' Meta_Data.deployment '.mat']),'ctd','-v7.3');
% save(fullfile(Meta_Data.CTDpath, ['alt_' Meta_Data.deployment '.mat']),'alt','-v7.3');

end


function [epsi,ctd,alt,act]=mod_som_read_epsi_raw(filename,Meta_Data)

c3515 = 42.914;

fid = fopen(filename);
fseek(fid,0,1);
frewind(fid);
str = fread(fid,'*char')';
fclose(fid);

% find pattern
% data format:
% EFEhex_timestamp,hex_length_block,hex_element_skipped,hex_voltage,hex_errorflagEFEDATA*chksum\r\n

%% Read EFE block depending on when file was created (block format has changed)

% NC added to switch how we read files depending on date they were created
% fileStruct = dir(filename);
% modifiedDate = fileStruct.datenum;
% NC 10/11/21 - Separated mod_som_read_epsi_files_v2.m and
% mod_som_read_epsi_files_v1.m. This one is v2;
efe_block_version = 'v2';


[ind_efe_start, ind_efe_end,ind_efe_tokens] = regexp(str,'\$EFE([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end','tokenExtents');

% if modifiedDate==datenum('31-Jan-2008 23:00:00')
%     efe_block_version = 'v2'; %raw file from SD card after 2021,4,1 that was not initialized
% elseif modifiedDate<datenum(2021,4,1)
%     efe_block_version = 'v1';
% elseif (modifiedDate>=datenum(2021,4,1) && modifiedDate<=datenum(2021,5,23))
%     efe_block_version = 'v2';
% elseif modifiedDate>=datenum(2021,5,23)
%     efe_block_version = 'v3';
% end
% 
% if (strcmp(str(ind_efe_start+30),"*")==0)
%     efe_block_version = 'v3';
% end

switch efe_block_version
    
    case 'v1'
        
        % token indices starts at ends at the square bracket
        % EFE[       1         ] [ 2  ]\r\n
        ind_efe_start=ind_efe_start+1;
        ind_efe_end=ind_efe_end+1;
        ind_efe_tokens=cellfun(@(x) (x+1),ind_efe_tokens,'un',0);
        % get the offsets to parse the str
        
        % define blocks header offset
        efe.header.strvalue  = 'EFE';
        efe.header.offset = 1;
        efe.header.length = length(efe.header.strvalue);
        
        efe.hextimestamp.strvalue  = "0000000000000000";
        efe.hextimestamp.length = strlength(efe.hextimestamp.strvalue);
        efe.hextimestamp.offset = efe.header.offset+efe.header.length;
        
        efe.hexlengthblock.strvalue  = "0000";
        efe.hexlengthblock.length = strlength(efe.hexlengthblock.strvalue);
        efe.hexlengthblock.offset = efe.hextimestamp.offset+ ...
            efe.hextimestamp.length+1; % +1 beacuse of the coma ","
        
        efe.hexelmntskip.strvalue  = "0000";
        efe.hexelmntskip.length = strlength(efe.hexelmntskip.strvalue);
        efe.hexelmntskip.offset = efe.hexlengthblock.offset+ ...
            efe.hexlengthblock.length+1; % +1 beacuse of the coma ","
        
        efe.hexerror.strvalue  = "0000";
        efe.hexerror.length = strlength(efe.hexerror.strvalue);
        efe.hexerror.offset = efe.hexelmntskip.offset+ ...
            efe.hexelmntskip.length+1; % +1 beacuse of the coma ","
        
        efe.chksum.strvalue = "FFFFF" ;
        efe.chksum.length   = strlength(efe.chksum.strvalue);
        
        % NC added to get laptop time (data before $EFE)
        efe.laptoptime.strvalue = '0000000000';
        efe.laptoptime.offset = -10;
        efe.laptoptime.length = 10;
        
        efe.data_offset = efe.hexerror.offset+efe.hexerror.length-1;
        
        % define some quantities
        
        efe.data.nchannels = 7;
        efe.data.sample_freq = 320;
        efe.data.sample_period = 1/efe.data.sample_freq;
        efe.data.bytes_per_channel = 3;
        efe.data.timestamp_length=8;
        efe.data.elementlength = efe.data.timestamp_length + ...
            efe.data.nchannels* efe.data.bytes_per_channel; % 8 bytes timestamps + 3 bytes ADC
        efe.data.nblocks = 160;
        
        efe.n_recs = numel(ind_efe_tokens);
        
        efe.data.length = efe.data.nblocks*efe.data.nchannels*efe.data.bytes_per_channel;
        
        efe.data.n_recs = numel(ind_efe_start);
        efe.time_sendout = NaN(efe.data.n_recs,1);
        efe.time_laptop = NaN(efe.data.n_recs,1);
        efe.bad_blocks = false(efe.data.n_recs,1);
        
    case 'v2'
        
        % token indices starts at ends at the square bracket
        % EFE[       1         ] [ 2  ]\r\n
        ind_efe_start=ind_efe_start+1;
        ind_efe_end=ind_efe_end+1;
        ind_efe_tokens=cellfun(@(x) (x+1),ind_efe_tokens,'un',0);
        
        % get the offsets to parse the str
        
        % define blocks header offset
        efe.header.strvalue  = 'EFE';
        efe.header.offset = 1;
        efe.header.length = length(efe.header.strvalue);
        
        efe.hextimestamp.strvalue  = "0000000000000000";
        efe.hextimestamp.length = strlength(efe.hextimestamp.strvalue);
        efe.hextimestamp.offset = efe.header.offset+efe.header.length;
        
        efe.hexlengthblock.strvalue  = "00000000";
        efe.hexlengthblock.length = strlength(efe.hexlengthblock.strvalue);
        efe.hexlengthblock.offset = efe.hextimestamp.offset+ ...
            efe.hextimestamp.length+1; % +1 beacuse of the coma ","
        
        efe.hexelmntskip.strvalue  = "00000000";
        efe.hexelmntskip.length = strlength(efe.hexelmntskip.strvalue);
        efe.hexelmntskip.offset = efe.hexlengthblock.offset+ ...
            efe.hexlengthblock.length+1; % +1 beacuse of the coma ","
        
        efe.hexerror.strvalue  = "0000";
        efe.hexerror.length = strlength(efe.hexerror.strvalue);
        efe.hexerror.offset = efe.hexelmntskip.offset+ ...
            efe.hexelmntskip.length+1; % +1 beacuse of the coma ","
        
        efe.chksum.strvalue = "FFFFF" ;
        efe.chksum.length   = strlength(efe.chksum.strvalue);
        
        % NC added to get laptop time (data before $EFE)
        efe.laptoptime.strvalue = '0000000000';
        efe.laptoptime.offset = -10;
        efe.laptoptime.length = 10;
        
        efe.data_offset = efe.hexerror.offset+efe.hexerror.length-1;
        
        % define some quantities
        
        efe.data.nchannels = 7;
        efe.data.sample_freq = 320;
        efe.data.sample_period = 1/efe.data.sample_freq;
        efe.data.bytes_per_channel = 3;
        efe.data.timestamp_length=8;
        efe.data.elementlength = efe.data.timestamp_length + ...
            efe.data.nchannels*efe.data.bytes_per_channel; % 8 bytes timestamps + 3 bytes ADC
        efe.data.nblocks = 160;
        
        efe.n_recs = numel(ind_efe_tokens);
        
        efe.data.length = efe.data.nblocks*efe.data.nchannels*efe.data.bytes_per_channel;
        
        efe.data.n_recs = numel(ind_efe_start);
        % NC added to keep track of laptop time
        efe.time_sendout = NaN(efe.data.n_recs,1);
        efe.time_laptop = NaN(efe.data.n_recs,1);
        efe.bad_blocks = false(efe.data.n_recs,1);
        
    case 'v3'
        
        % token indices starts at ends at the square bracket
        % EFE[       1         ] [ 2  ]\r\n
        ind_efe_tokens=cellfun(@(x) (x+1),ind_efe_tokens,'un',0);
        
        % get the offsets to parse the str
        efe.sync.strvalue='$';
        efe.sync.offset=0;
        efe.sync.length=1;
        % define blocks header offset
        efe.header.strvalue  = 'EFE4';
        efe.header.length = length(efe.header.strvalue)-1;
        efe.header.offset = efe.sync.offset+efe.sync.length+1;
        
        efe.hextimestamp.strvalue  = "0000000000000000";
        efe.hextimestamp.length = strlength(efe.hextimestamp.strvalue)-1;
        efe.hextimestamp.offset = efe.header.offset+efe.header.length+1;
        
        efe.hexlengthblock.strvalue  = "00000000";
        efe.hexlengthblock.length = strlength(efe.hexlengthblock.strvalue)-1;
        efe.hexlengthblock.offset = efe.hextimestamp.offset+ ...
            efe.hextimestamp.length+1;
        
        efe.headerchecksum.strvalue = "*FF";
        efe.headerchecksum.length   = strlength(efe.headerchecksum.strvalue)-1;
        efe.headerchecksum.offset   = efe.hexlengthblock.offset+ ...
            efe.hexlengthblock.length+1;
        
        efe.chksum.strvalue = "FFFFF" ;
        efe.chksum.length   = strlength(efe.chksum.strvalue);
        
        % NC added to get laptop time (data before $EFE)
        efe.laptoptime.strvalue = '0000000000';
        efe.laptoptime.offset = -10;
        efe.laptoptime.length = 10;
        
        
        efe.data_offset = efe.headerchecksum.offset + ...
            efe.headerchecksum.length;
        
        % define some quantities
        
        efe.data.nchannels = 7;
        efe.data.sample_freq = 320;
        efe.data.sample_period = 1/efe.data.sample_freq;
        efe.data.bytes_per_channel = 3;
        efe.data.timestamp_length=8;
        efe.data.elementlength = efe.data.timestamp_length + ...
            efe.data.nchannels*efe.data.bytes_per_channel; % 8 bytes timestamps + 3 bytes ADC
        efe.data.nblocks = 160;
        
        efe.n_recs = numel(ind_efe_tokens);
        
        efe.data.length = efe.data.nblocks*efe.data.nchannels*efe.data.bytes_per_channel;
        
        efe.data.n_recs = numel(ind_efe_start);
        % NC added to keep track of laptop time
        efe.time_sendout = NaN(efe.data.n_recs,1);
        efe.time_laptop = NaN(efe.data.n_recs,1);
        efe.bad_blocks = false(efe.data.n_recs,1);
        
end

switch efe_block_version
    case {'v1','v2'}
        get_inds = @(x)(x.offset:x.offset+x.length-1);
        sparse_header=@(header,x) (hex2dec(header(x.offset: ...
            x.offset+ ...
            x.length-1)));
    case {'v3'}
        get_inds = @(x)(x.offset:x.offset+x.length);
        get_inds_alti = @(x)(x.offset:x.offset+x.length-1);
        sparse_header=@(header,x) (hex2dec(header(x.offset: ...
            x.offset+ ...
            x.length)));
end

%% read the EFE record inside str
if isempty(ind_efe_start)
    disp('no epsi data')
    epsi=[];
else
    
    for i = 1:efe.data.n_recs
        
        %fprintf("nb block %i \n",i);
        
        % get the header of the block
        header= str(ind_efe_start(i)+(0:efe.data_offset-1));
        switch efe_block_version
            case {'v1','v2'}
                % parse the header of the block
                efe.hextimestamp.value=hex2dec(header(get_inds(efe.hextimestamp)));
                efe.hexlengthblock.value=hex2dec(header(get_inds(efe.hexlengthblock)));
                efe.hexelmntskip.value=hex2dec(header(get_inds(header,efe.hexelmntskip)));
                efe.hexerror.value=hex2dec(header(get_inds(header,efe.hexerror)));
                
                % NC added to get laptop time
                % Get laptop time string (10 characters before $EFE header)
                laptoptime_str = str(ind_efe_start(i)-1+...
                    (efe.laptoptime.offset:...
                    efe.laptoptime.offset+efe.laptoptime.length-1));
                efe.laptoptime.value = str2double(laptoptime_str); %laptop time in hundredths of seconds since ???
                
                % compare the length of block from regexp and from the header
                Lregexp=ind_efe_end(i)-ind_efe_start(i);
            case {'v3'}
                % parse the header of the block
                efe.hextimestamp.strvalue=header(efe.hextimestamp.offset+ ...
                    (0:efe.hextimestamp.length));
                efe.hextimestamp.value=hex2dec(header(get_inds(efe.hextimestamp)));
                
                efe.hexlengthblock.strvalue=header(efe.hexlengthblock.offset+ ...
                    (0:efe.hexlengthblock.length));
                efe.hexlengthblock.value=hex2dec(header(get_inds(efe.hexlengthblock)));
                
                % compare the length of block from regexp and from the header
                Lregexp=ind_efe_end(i)-ind_efe_start(i)-efe.data_offset-5;
        end
        
        if(efe.hexlengthblock.value~=Lregexp+1)
            % NC added to track bad blocks
            efe.bad_blocks(i) = true;
            fprintf("block %i: bad block\n",i)
        else
            
            switch efe_block_version
                case {'v1','v2'}
                    Lefeblock=efe.hexlengthblock.value- ...
                        efe.data_offset- ...
                        efe.chksum.length-1;
                case {'v3'}
                    Lefeblock=efe.hexlengthblock.value;
            end
            
            nb_element=(Lefeblock)./efe.data.elementlength;
            
            raw_bytes=uint32(str(ind_efe_start(i)+efe.data_offset+ ...
                (0:Lefeblock-1)));
            efe.data.raw_bytes{i}=reshape(raw_bytes,efe.data.elementlength,nb_element).';        
            
            switch efe_block_version
                case {'v1','v2'}
                    % NC changed definition for efe.checksum.data to work for v3, but
                    % it might break it for v1 and v2. Stop here and make sure the
                    % input to hex2dec is two characters. If not, go up to where chksum
                    % offset is defined and change it (and anything else that uses
                    % chksum offset that might need to be changed
                    str(ind_efe_start(i)+ ...
                        efe.data_offset+Lefeblock+(1:2))
                    pause
                case 'v3'
                    
                    efe.checksum.data(i) = hex2dec(str(ind_efe_start(i)+ ...
                        efe.data_offset+Lefeblock+(1:2)));
            end
            
            % NC added to track hextimestamp (time_sendout) and laptop time
            % (time_laptop)
            efe.time_sendout(i) = efe.hextimestamp.value;
            try
                efe.time_laptop(i) = efe.laptoptime.value./100;
            catch
                efe.time_laptop(i) = nan;
            end
        end
        
    end
    
    
    %% clean and concat efe.data.raw_bytes
    idnull=cellfun(@isempty,efe.data.raw_bytes);
    efe.data.raw_bytes=cell2mat(efe.data.raw_bytes(~idnull).');
    
    %% get epsi timestamp in the epsi structure
    epsi.timestamp=zeros(size(efe.data.raw_bytes,1),1);
    powers = [0:7].';
    
    epsi.timestamp  = double(uint64(efe.data.raw_bytes(:,1:8)))*256.^powers; %NC - timestamp in milliseconds since 1970 is too big for uint32
   
    %% get ADC data
    efe.data.raw = efe.data.raw_bytes(:,efe.data.timestamp_length+1:3:end)*256^2+ ...
        efe.data.raw_bytes(:,efe.data.timestamp_length+2:3:end)*256+ ...
        efe.data.raw_bytes(:,efe.data.timestamp_length+3:3:end);
    
       
    %% convert count 2 volt
    
    bit_counts = 24;
    gain = 1;
    acc_offset = 1.8/2;
    acc_factor = 0.5;
    channels=Meta_Data.PROCESS.channels;
    nb_channels=length(channels);
    
    for cha=1:nb_channels
        wh_channel=channels{cha};
        epsi.([wh_channel '_count']) = efe.data.raw(:,cha);
    end
    
    Unipolar=@(FR,data) (FR/gain*double(data)/2.^(bit_counts));
    Bipolar=@(FR,data) (FR/gain*(double(data)/2.^(bit_counts-1)-1));
    for cha=1:nb_channels
        wh_channel=channels{cha};
        FR=Meta_Data.AFE.(wh_channel).full_range;
        switch Meta_Data.AFE.(wh_channel).ADCconf
            case {'Bipolar','bipolar'}
                epsi.([wh_channel '_volt'])=Bipolar(FR,epsi.([wh_channel '_count']));
            case {'Unipolar','unipolar'}
                epsi.([wh_channel '_volt'])=Unipolar(FR,epsi.([wh_channel '_count']));
                
        end
        
        switch Meta_Data.AFE.(wh_channel).full_range
            case 1.8 %case for acceleration channels
                epsi.([wh_channel '_g']) = (epsi.([wh_channel '_volt'])-acc_offset)/acc_factor;
                epsi = rmfield(epsi,[wh_channel '_volt']);
        end
    end
    
    epsi_fields = fieldnames(epsi);
    for i  = 1:numel(epsi_fields)
        epsi.(epsi_fields{i}) = reshape(epsi.(epsi_fields{i}),[],1);
    end
    
    % If timestamp has values like 1.6e12, it is in milliseconds since Jan
    % 1, 1970. Otherwise it's in milliseconds since the start of the record
    if epsi.timestamp(i)>1e9
        seconds1970 = epsi.timestamp./1000;
        days1970 = seconds1970/(24 * 60 * 60);
        offset1970_days = datenum(1970, 1, 1) - datenum(0000, 1, 0);
        offset1970_seconds = seconds(days(offset1970_days));
        epsi.epsitime = seconds1970 + offset1970_seconds;
        epsi.epsidnum = days1970 + offset1970_days;
    else
        epsi.epsitime = epsi.timestamp/1000;
    end
    
    % NC add the efe timestamp fields
    % There is one efe time sample for every 80 epsi samples ($EFE blocks
    % contain 80 data samples. Replicate the efe timestamps to match size of
    % data
    fieldsToExpand = {'bad_blocks','time_sendout','time_laptop'};
    for iField=1:length(fieldsToExpand)
        a = repmat(efe.(fieldsToExpand{iField})(:).',nb_element,1);
        b = reshape(a,efe.n_recs*nb_element,1);
        epsi.(['efe_' fieldsToExpand{iField}]) = b;
    end
    block_start = repmat([1;nan(nb_element-1,1)],1,efe.n_recs);
    epsi.efe_block_start = reshape(block_start,efe.n_recs*nb_element,1);
    
end


%% SBE 49
switch efe_block_version
    case {'v1','v2'}
        [ind_sbe_start,ind_sbe_stop, ind_sbe_tokens] = regexp(str,'\$S49([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end','tokenExtents');
    case 'v3'
        [ind_sbe_start,ind_sbe_stop, ind_sbe_tokens] = regexp(str,'\$SB49([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end','tokenExtents');
end

if isempty(ind_sbe_start)
    [ind_sbe_start,ind_sbe_stop, ind_sbe_tokens] = regexp(str,'\$SB41([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end','tokenExtents');
    if isempty(ind_sbe_start)
        disp('no ctd data')
        ctd=[];
    end
else
    
    if (strcmp(str(ind_sbe_start+30),"*")==0)
        sbe_block_version = 'v3';
    else
        sbe_block_version = 'v2';
    end
    
    sbe.data.n_block  = numel(ind_sbe_start);
    sbe.data.n_recs   = numel(ind_sbe_start)*Meta_Data.CTD.sample_per_record;
    switch sbe_block_version
        case {'v1','v2'}
            
            sbe.header.strvalue    = 'S49'; % or SBE41
            sbe.header.offset = 2;
            sbe.header.length = strlength(sbe.header.strvalue);
            
            sbe.hextimestamp.strvalue  = "0000000000000000";
            sbe.hextimestamp.length    = strlength(sbe.hextimestamp.strvalue);
            sbe.hextimestamp.offset    = sbe.header.offset+sbe.header.length;
            
            sbe.hexlengthblock.strvalue  = "0000";
            sbe.hexlengthblock.length    = strlength(sbe.hexlengthblock.strvalue);
            sbe.hexlengthblock.offset    = sbe.hextimestamp.offset+ ...
                sbe.hextimestamp.length+1; % +1 beacuse of the coma ","
            
            sbe.hexelmntskip.strvalue  = "0000";
            sbe.hexelmntskip.length    = strlength(sbe.hexelmntskip.strvalue);
            sbe.hexelmntskip.offset    = sbe.hexlengthblock.offset+ ...
                sbe.hexlengthblock.length+1; % +1 beacuse of the coma ","
            
            
            sbe.data_offset = sbe.hexelmntskip.offset+sbe.hexelmntskip.length-1;
        case 'v3'
            
            ind_sbe_tokens=cellfun(@(x) (x+1),ind_sbe_tokens,'un',0);
            
            % get the offsets to parse the str
            sbe.sync.strvalue='$';
            sbe.sync.offset=0;
            sbe.sync.length=1;
            % define blocks header offset
            sbe.header.strvalue  = 'SB49';
            sbe.header.length = length(sbe.header.strvalue)-1;
            sbe.header.offset = sbe.sync.offset+sbe.sync.length+1;
            
            sbe.hextimestamp.strvalue  = "0000000000000000";
            sbe.hextimestamp.length = strlength(sbe.hextimestamp.strvalue)-1;
            sbe.hextimestamp.offset = sbe.header.offset+sbe.header.length+1;
            
            sbe.hexlengthblock.strvalue  = "00000000";
            sbe.hexlengthblock.length = strlength(sbe.hexlengthblock.strvalue)-1;
            sbe.hexlengthblock.offset = sbe.hextimestamp.offset+ ...
                sbe.hextimestamp.length+1;
            
            sbe.headerchecksum.strvalue = "*FF";
            sbe.headerchecksum.length   = strlength(sbe.headerchecksum.strvalue)-1;
            sbe.headerchecksum.offset   = sbe.hexlengthblock.offset+ ...
                sbe.hexlengthblock.length+1;
            
            sbe.chksum.strvalue = "FFFFF" ;
            sbe.chksum.length   = strlength(sbe.chksum.strvalue);
            
            % NC added to get laptop time (data before $EFE)
            sbe.laptoptime.strvalue = '0000000000';
            sbe.laptoptime.offset = -10;
            sbe.laptoptime.length = 10;
            
            
            sbe.data_offset = sbe.headerchecksum.offset + ...
                sbe.headerchecksum.length;
            
    end
    % TODO fix the SBE name bug
    switch Meta_Data.CTD.name
        case{"SBE49","SBE","S49","SB49"}
            sbe.data.format = 'eng';
            sbe.data.length=22;
            sbe.data.sample_freq = 16;
            sbe.cal        = Meta_Data.CTD.cal;
            ctd.P_raw      = NaN(sbe.data.n_recs,1);
            ctd.T_raw      = NaN(sbe.data.n_recs,1);
            ctd.S_raw      = NaN(sbe.data.n_recs,1);
            ctd.C_raw      = NaN(sbe.data.n_recs,1);
            ctd.PT_raw     = NaN(sbe.data.n_recs,1);
            
        case{"SBE41","S41","SB41"}
            sbe.data.format = 'PTS';
            sbe.data.length=28;
            sbe.data.sample_freq = 1;
    end
    
    sbe.data.timestamp_length=16;
    sbe.data.sample_period  = 1/sbe.data.sample_freq;
    sbe.block_time          = NaN(sbe.data.n_block,1);
    sbe.checksum            = uint8(zeros(sbe.data.n_recs,1));
    
    ctd.timestamp = zeros(sbe.data.n_recs,1);
    ctd.P         = NaN(sbe.data.n_recs,1);
    ctd.T         = NaN(sbe.data.n_recs,1);
    ctd.S         = NaN(sbe.data.n_recs,1);
    ctd.C         = NaN(sbe.data.n_recs,1);
    
    n_rec=0;
    for i = 1:sbe.data.n_block
        % might want to check the length of data too!!! right now I am skipping
        % that step
        tmp_sbe_block = str(ind_sbe_start(i):ind_sbe_stop(i)); % +2 becasue now the header is SBE41 or SBE49 and not SBE
        switch sbe_block_version
            case {'v1','v2'}
                sbe.hextimestamp.value=hex2dec(tmp_sbe_block(get_inds(sbe.hextimestamp)));
                sbe.hexlengthblock.value=hex2dec(tmp_sbe_block(get_inds(sbe.hexlengthblock)));
                sbe.hexelmntskip.value=hex2dec(tmp_sbe_block(get_inds(sbe.hexelmntskip)));
                
                sbe.block_time(i)  = sbe.hextimestamp.value./1000;
            case 'v3'
                sbe.hextimestamp.value=hex2dec(tmp_sbe_block(get_inds(sbe.hextimestamp)));
                sbe.hexlengthblock.value=hex2dec(tmp_sbe_block(get_inds(sbe.hexlengthblock)));
                
                %                 % compare the length of block from regexp and from the header
                %                 Lregexp=ind_sbe_end(i)-ind_sbe_start(i)-sbe.data_offset-5;
                
        end
        tmp_data_ctd=tmp_sbe_block(sbe.data_offset+1:end-5);
        
        if (length(tmp_data_ctd)~=(24+16)*Meta_Data.CTD.sample_per_record)
            fprintf("not enough SBE sample in block %i\r\n",i)
        else
            %
            tmp_data_ctd=reshape(tmp_data_ctd,16+24,Meta_Data.CTD.sample_per_record).';
            
            for j=1:Meta_Data.CTD.sample_per_record
                n_rec=n_rec+1;
                element_ctd=tmp_data_ctd(j,1:end-2);
                ctd.timestamp(n_rec) = hex2dec(element_ctd(1:16));
                rec_ctd=element_ctd(sbe.data.timestamp_length+1:end);
                %             rec_ctd=rec_ctd(rec_ctd~=' ');
                
                switch sbe.data.format
                    case 'PTS'
                        data_split   = strsplit(rec_ctd,',');
                        ctd.P(n_rec)        = str2double(data_split{1});
                        ctd.T(n_rec)        = str2double(data_split{2});
                        ctd.S(n_rec)        = str2double(data_split{3});
                        ctd.C(n_rec)        = NaN;
                    case 'eng'
                        raw_sample   = rec_ctd;
                        ctd.T_raw(n_rec) = hex2dec(raw_sample(:,1:6));
                        ctd.C_raw(n_rec) = hex2dec(raw_sample(:,(1:6)+6));
                        ctd.P_raw(n_rec) = hex2dec(raw_sample(:,(1:6)+12));
                        ctd.PT_raw(n_rec) = hex2dec(raw_sample(:,(1:4)+18));
                end
            end
        end
    end
    
    % If timestamp has values like 1.6e12, it is in milliseconds since Jan
    % 1, 1970. Otherwise it's in milliseconds since the start of the record
    if ctd.timestamp(n_rec)>1e9
        seconds1970 = ctd.timestamp./1000;
        days1970 = seconds1970/(24 * 60 * 60);
        offset1970_days = datenum(1970, 1, 1) - datenum(0000, 1, 0);
        offset1970_seconds = seconds(days(offset1970_days));
        ctd.ctdtime = seconds1970 + offset1970_seconds;
        ctd.ctddnum = days1970 + offset1970_days;
    else
        ctd.ctdtime = ctd.timestamp./1000;
    end
    
    %%
    switch sbe.data.format
        case 'eng'
            ctd = sbe49_ascii_get_temperature(ctd,sbe);
            ctd = sbe49_ascii_get_pressure(ctd,sbe);
            ctd = sbe49_ascii_get_conductivity(ctd,sbe);
            
            ctd.S    = real(sw_salt(ctd.C*10./c3515,ctd.T,ctd.P)); %Every once in a while, S comes out imaginary, I think only when SBE is on deck. 
            ctd.sig  = sw_pden(ctd.S,ctd.T,ctd.P,0);
            ctd.dPdt = [0; diff(ctd.P)./diff(ctd.ctdtime)];
            
            
            
    end
    
    
end


%% ALTI

[ind_alti_start,ind_alti_stop, ind_alti_tokens] = regexp(str,'\$ALT([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end','tokenExtents');
if isempty(ind_alti_start)
    disp('no alti data')
    alt=[];
else
    
    if (strcmp(str(ind_alti_start+30),"*")==0)
        alti_block_version = 'v3';
    else
        alti_block_version = 'v1';
    end
    
    alti.data.n_block  = numel(ind_alti_start);
    alti.data.n_recs   = numel(ind_alti_stop);
    
    alti.soundSpeed    = 1500;
    
    switch alti_block_version
        case {"v1","v2"}
            alti.header.strvalue    = 'ALT'; % or SBE41
            alti.header.offset = 2;
            alti.header.length = strlength(alti.header.strvalue);
            
            alti.hextimestamp.strvalue  = "0000000000000000";
            alti.hextimestamp.length    = strlength(alti.hextimestamp.strvalue);
            alti.hextimestamp.offset    = alti.header.offset+alti.header.length;
            
            alti.data_offset = alti.hextimestamp.offset+alti.hextimestamp.length+1;
            
        case 'v3'
            
            % get the offsets to parse the str
            alti.sync.strvalue='$';
            alti.sync.length=strlength(alti.sync.strvalue);
            alti.sync.offset=1;
            
            % define blocks header offset
            alti.header.strvalue  = 'ALTI';
            alti.header.length = strlength(alti.header.strvalue);
            alti.header.offset = alti.sync.offset+alti.sync.length;
            
            alti.hextimestamp.strvalue  = "0000000000000000";
            alti.hextimestamp.length = strlength(alti.hextimestamp.strvalue);
            alti.hextimestamp.offset = alti.header.offset+alti.header.length;
            
            alti.hexlengthblock.strvalue  = "00000000";
            alti.hexlengthblock.length = strlength(alti.hexlengthblock.strvalue);
            alti.hexlengthblock.offset = alti.hextimestamp.offset+ ...
                alti.hextimestamp.length;
            
            alti.headerchecksum.strvalue = "*FF";
            alti.headerchecksum.length   = strlength(alti.headerchecksum.strvalue);
            alti.headerchecksum.offset   = alti.hexlengthblock.offset+ ...
                alti.hexlengthblock.length;
            
            alti.chksum.strvalue = "FFFFF" ;
            alti.chksum.length   = strlength(alti.chksum.strvalue);
            
            % NC added to get laptop time (data before $EFE)
            alti.laptoptime.strvalue = '0000000000';
            alti.laptoptime.offset = -10;
            alti.laptoptime.length = 10;
            
            alti.data_offset = alti.hextimestamp.offset + ...
                alti.hextimestamp.length;
            
    end
    
    alt.alttime = [];
    alt.dst = [];
    for i = 1:alti.data.n_block
        % might want to check the length of data too!!! right now I am skipping
        % that step
        tmp_alt_block = str(ind_alti_start(i):ind_alti_stop(i)); % +2 becasue now the header is SBE41 or SBE49 and not SBE
        %alti.hextimestamp.value=sparse_header(tmp_alt_block,alti.hextimestamp);
        alti.hextimestamp.value = hex2dec(tmp_alt_block(get_inds_alti(alti.hextimestamp)));
        alt.timestamp(i,1)  = alti.hextimestamp.value;
        tmp_data_alt=tmp_alt_block(alti.data_offset:end-5);
        
        data_split = strsplit(tmp_data_alt,'*');
        alt.dst(i,1)=str2double(data_split{1})*1e-5*alti.soundSpeed;
        
    end
    
    % If timestamp has values like 1.6e12, it is in milliseconds since Jan
    % 1, 1970. Otherwise it's in milliseconds since the start of the record
    if alt.timestamp(i)>1e9
        seconds1970 = alt.timestamp./1000;
        days1970 = seconds1970/(24 * 60 * 60);
        offset1970_days = datenum(1970, 1, 1) - datenum(0000, 1, 0);
        offset1970_seconds = seconds(days(offset1970_days));
        alt.alttime = seconds1970 + offset1970_seconds;
        alt.altdnum = days1970 + offset1970_days;
    else
        alt.alttime = alt.timestamp/1000;
    end
        
end

%% ACTU

[ind_actu_start,ind_actu_stop, ind_actu_tokens] = regexp(str,'\$ACTU([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end','tokenExtents');
if isempty(ind_actu_start)
    disp('no actu data')
    act=[];
else
    disp('actu data')
    
    actu.data.n_block  = numel(ind_actu_start);
    actu.data.n_recs   = numel(ind_actu_stop);
    
    actu.header.strvalue    = 'ACT'; % or SBE41
    actu.header.offset = 2;
    actu.header.length = strlength(actu.header.strvalue);
    
    actu.hextimestamp.strvalue  = "0000000000000000";
    actu.hextimestamp.length    = strlength(actu.hextimestamp.strvalue);
    actu.hextimestamp.offset    = actu.header.offset+actu.header.length;
    
    actu.data_offset = actu.hextimestamp.offset+actu.hextimestamp.length+1;
    
    act.alttime = [];
    act.dst = [];
    for i = 1:actu.data.n_block
        % might want to check the length of data too!!! right now I am skipping
        % that step
        tmp_act_block = str(ind_actu_start(i):ind_actu_stop(i)); % +2 becasue now the header is SBE41 or SBE49 and not SBE
        actu.hextimestamp.value=hex2dec(tmp_act_block(get_inds(actu.hextimestamp)));
        
        act.timestamp(i,1)  = actu.hextimestamp.value;
        tmp_data_act=tmp_act_block(actu.data_offset:end-5);
        act.dst(i,1)=str2double(tmp_data_act)*1e-5*alti.soundSpeed;
    end
    
    % If timestamp has values like 1.6e12, it is in milliseconds since Jan
    % 1, 1970. Otherwise it's in milliseconds since the start of the record
    if act.timestamp(i)>1e9
        seconds1970 = act.timestamp./1000;
        days1970 = seconds1970/(24 * 60 * 60);
        offset1970 = datenum(1970, 1, 1) - datenum(0000, 1, 0);
        act.altdnum = days1970 + offset1970;
    else
        act.acttime = act.timestamp/1000;
    end
    
end


end

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

% check for empty efe data
old_varargin=varargin;
old_nargin=nargin;
empty_efe=cellfun(@isempty,varargin);
varargin=varargin(~empty_efe);
new_nargin=numel(varargin);

%sort the files
start_time_file=cellfun(@(x) x.epsitime(1), varargin);
[~,I]=sort(start_time_file);
varargin=varargin(I);

epsi_fields = fieldnames(varargin{1});
for i = 2:new_nargin
    tmp_fields = fieldnames(varargin{i});
    for j = 1:numel(tmp_fields)
        if ~ismember(tmp_fields{j},epsi_fields)
            epsi_fields{end+1} = tmp_fields{j};
        end
    end
end

% epsi_sub_fields = cell(size(epsi_fields));
%
% for i = 1:numel(epsi_fields)
%     epsi_sub_fields{i} = fieldnames(varargin{1}.(epsi_fields{i}));
%
%     for j = 2:nargin
%         tmp_fields = fieldnames(varargin{j}.(epsi_fields{i}));
%         for k = 1:numel(tmp_fields)
%             if ~ismember(tmp_fields{k},epsi_sub_fields{i})
%                 epsi_sub_fields{i}{end+1} = tmp_fields{k};
%             end
%         end
%     end
% end

for i=1:(length(epsi_fields))
    %header field
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
end
end


function ctd = mod_combine_ctd(varargin)
% mod_combine_epsi - combines epsi data files in MATLAB format that was
% converted using mod_read_epsi_raw
%
% mod_combine_epsi(epsi1,epsi2,epsi3,...) returns a EPSI structure of variables described
% for MET data files
%
% Written 2018/10/15 - San Nguyen stn 004@ucsd.edu

if nargin < 1
    ctd = [];
    return;
end

if nargin == 1
    ctd = varargin{1};
    if length(ctd) == 1
        return;
    end
    evalstr = 'ctd = mod_combine_ctd(';
    for i=1:(length(ctd)-1)
        evalstr = strcat(evalstr, 'ctd(', num2str(i), '),');
    end
    evalstr = strcat(evalstr, 'ctd(', num2str(length(ctd)), '));');
    eval(evalstr);
    return
end

% check emtpy structures
old_varargin=varargin;
old_nargin=nargin;
empty_ctd=cellfun(@isempty,varargin);
varargin=varargin(~empty_ctd);
new_nargin=numel(varargin);
start_time_file=cellfun(@(x) x.ctdtime(1), varargin);
[~,I]=sort(start_time_file);
varargin=varargin(I);

% grab the data field name
ctd_fields = fieldnames(varargin{1});

for i=1:length(ctd_fields)
    evalstr = strcat('ctd.', ctd_fields{i}, '= [');
    for j=1:new_nargin
        if ~isfield(varargin{j},(ctd_fields{i}))
            varargin{j}.(ctd_fields{i}) = NaN(size(varargin{j}.ctdtime));
        end
        evalstr = strcat(evalstr, 'varargin{', num2str(j), '}.', ctd_fields{i},';');
    end
    evalstr = strcat(evalstr(1:end-1), '];');
    eval(evalstr);
end

end

function alt = mod_combine_alt(varargin)
% mod_combine_epsi - combines epsi data files in MATLAB format that was
% converted using mod_read_epsi_raw
%
% mod_combine_epsi(epsi1,epsi2,epsi3,...) returns a EPSI structure of variables described
% for MET data files
%
% Written 2018/10/15 - San Nguyen stn 004@ucsd.edu

if nargin < 1
    alt = [];
    return;
end

if nargin == 1
    alt = varargin{1};
    if length(alt) == 1
        return;
    end
    evalstr = 'alt = mod_combine_alt(';
    for i=1:(length(alt)-1)
        evalstr = strcat(evalstr, 'alt(', num2str(i), '),');
    end
    evalstr = strcat(evalstr, 'alt(', num2str(length(alt)), '));');
    eval(evalstr);
    return
end

% check emtpy structures
old_varargin=varargin;
old_nargin=nargin;
empty_alt=cellfun(@isempty,varargin);
varargin=varargin(~empty_alt);
new_nargin=numel(varargin);
start_time_file=cellfun(@(x) x.alttime(1), varargin);
[~,I]=sort(start_time_file);
varargin=varargin(I);

% grab the data field name
alt_fields = fieldnames(varargin{1});

for i=1:length(alt_fields)
    evalstr = strcat('alt.', alt_fields{i}, '= [');
    for j=1:new_nargin
        if ~isfield(varargin{j},(alt_fields{i}))
            varargin{j}.(alt_fields{i}) = NaN(size(varargin{j}.alttime));
        end
        evalstr = strcat(evalstr, 'varargin{', num2str(j), '}.', alt_fields{i},';');
    end
    evalstr = strcat(evalstr(1:end-1), '];');
    eval(evalstr);
end

end




%%
%  reads and apply calibration to the temperature data
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
