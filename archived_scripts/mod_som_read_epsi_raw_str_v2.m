function epsi=mod_som_read_epsi_raw_str_v2(str,Meta_Data)
%%

% find pattern
% data format:
% EFEhex_timestamp,hex_length_block,hex_element_skipped,hex_voltage,hex_errorflagEFEDATA*chksum\r\n

% token indices starts at ends at the square bracket
% EFE[       1         ] [ 2  ]\r\n
[ind_efe_start, ind_efe_end,ind_efe_tokens] = regexp(str,'EFE([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end','tokenExtents');

%% get the offsets to parse the str

% define blocks header offset
efe.header.strvalue  = 'EFE';
efe.header.offset = 0;
efe.header.length = length(efe.header.strvalue);

efe.hextimestamp.strvalue  = "0000000000000000";
efe.hextimestamp.length = strlength(efe.hextimestamp.strvalue);
efe.hextimestamp.offset = efe.header.length+1; 

efe.hexlengthblock.strvalue  = "0000";
efe.hexlengthblock.length = strlength(efe.hexlengthblock.strvalue);
efe.hexlengthblock.offset = efe.hextimestamp.offset+ ...
                            efe.hextimestamp.length+1; % +1 beacuse of the coma "," 

efe.hexelmntskip.strvalue  = "0000";
efe.hexelmntskip.length = strlength(efe.hexelmntskip.strvalue);
efe.hexelmntskip.offset = efe.hexlengthblock.offset+ ...
                          efe.hexlengthblock.length+1; % +1 beacuse of the coma "," 

efe.hexvoltage.strvalue  = "0000";
efe.hexvoltage.length = strlength(efe.hexvoltage.strvalue);
efe.hexvoltage.offset = efe.hexelmntskip.offset+ ...
                        efe.hexelmntskip.length+1; % +1 beacuse of the coma "," 

efe.hexerror.strvalue  = "0000";
efe.hexerror.length = strlength(efe.hexerror.strvalue);
efe.hexerror.offset = efe.hexvoltage.offset+ ...
                      efe.hexvoltage.length+1; % +1 beacuse of the coma "," 



efe.chksum.strvalue = "*FF\r\n" ;
efe.chksum.length   = strlength(efe.chksum.strvalue); 


efe.data_offset = efe.hexerror.offset+efe.hexerror.length-1;

%% define some quantities

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
efe.time = NaN(efe.data.n_recs,1);
efe.status = NaN(efe.data.n_recs,1);

sparse_header=@(header,x) (hex2dec(header(x.offset: ...
                                   x.offset+ ...
                                   x.length-1)));
%% read the EFE record inside str
for i = 1:efe.data.n_recs
    % get the header of the block
    header= str(ind_efe_start(i):ind_efe_start(i)+efe.data_offset-1);
    % parse the header of the block
    efe.hextimestamp.value=sparse_header(header,efe.hextimestamp);
    efe.hexlengthblock.value=sparse_header(header,efe.hexlengthblock);
    efe.hexelmntskip.value=sparse_header(header,efe.hexelmntskip);
    efe.hexvoltage.value=sparse_header(header,efe.hexvoltage);
    efe.hexerror.value=sparse_header(header,efe.hexerror);
    
    % compare the length of block from regexp and from the header
    Lregexp=ind_efe_end(i)-ind_efe_start(i);
    if(efe.hexlengthblock.value~=Lregexp+1)
        fprintf("block %i: bad block\n",i)
    else
        Lefeblock=efe.hexlengthblock.value-efe.data_offset+1-efe.chksum.length+1;
        
        nb_element=Lefeblock./efe.data.elementlength;
        
        %     tmp_data_u8 = uint8(str(ind_efe_start(i)+(0:Lregexp-efe.chksum.length+2)));
        
        raw_bytes=uint32(str(ind_efe_start(i)+efe.data_offset+ ...
            (0:Lefeblock-1)));
        efe.data.raw_bytes{i}=reshape(raw_bytes,efe.data.elementlength,nb_element).';
        
        
        efe.checksum.data(i) = hex2dec(str(ind_efe_start(i)+ ...
            efe.data_offset+Lefeblock+1+(0:1)));
        
        %     for j = 1:numel(tmp_data_u8) %% check sum is not working yet
        %         efe.checksum.calculate = bitxor(efe.checksum.calculate,tmp_data_u8(j));
        %     end
        efe.time(i) = efe.hextimestamp.value;
    end
    
end

%% clean and concat efe.data.raw_bytes
idnull=cellfun(@isempty,efe.data.raw_bytes);
efe.data.raw_bytes=cell2mat(efe.data.raw_bytes(~idnull).');

%% get epsi timestamp in the epsi structure
epsi.timestamp=zeros(size(efe.data.raw_bytes,1),1);
for i=efe.data.timestamp_length:-1:1
    epsi.timestamp  =  epsi.timestamp + ...
                                         double(efe.data.raw_bytes(:,i)*256^(i-1)); 
end

%% get ADC data 
efe.data.raw = efe.data.raw_bytes(:,efe.data.timestamp_length+1:3:end)*256^2+ ...
               efe.data.raw_bytes(:,efe.data.timestamp_length+2:3:end)*256+ ...
               efe.data.raw_bytes(:,efe.data.timestamp_length+3:3:end);

%% convert count 2 volt
% acc_offset = 1.8/2;
% acc_factor = 0.5;
% Meta_Data.epsi.t1.full_range=2.5;
% Meta_Data.epsi.t2.full_range=2.5;
% Meta_Data.epsi.s1.full_range=2.5;
% Meta_Data.epsi.s2.full_range=2.5;
% Meta_Data.epsi.a1.full_range=2.5;
% Meta_Data.epsi.a1.full_range=1.8;
% Meta_Data.epsi.a2.full_range=1.8;
% Meta_Data.epsi.a3.full_range=1.8;

bit_counts = 24;
gain = 1;
channels={'t1','t2','s1','s2','a1','a2','a3'};
nb_channels=length(channels);

for cha=1:nb_channels
    wh_channel=channels{cha};
    epsi.([wh_channel '_count']) = efe.data.raw(:,cha);
end

Unipolar=@(FR,data) (FR/gain*double(data)/2.^(bit_counts)); 
Bipolar=@(FR,data) (FR/gain*(double(data)/2.^(bit_counts-1)-1)); 
for cha=1:nb_channels
    wh_channel=channels{cha};
    FR=Meta_Data.epsi.(wh_channel).full_range;
    if ~strcmp(wh_channel,'c')
        switch Meta_Data.epsi.(wh_channel).ADCconf
            case {'Bipolar','bipolar'}
                epsi.(wh_channel)=Bipolar(FR,epsi.([wh_channel '_count']));
            case {'Unipolar','unipolar'}
                epsi.(wh_channel)=Unipolar(FR,epsi.([wh_channel '_count']));
                
        end
        
%         switch wh_channel
%             case 'a1'
%                 epsi.a1 = (epsi.a1-acc_offset)/acc_factor;
%             case 'a2'
%                 epsi.a2 = (epsi.a2-acc_offset)/acc_factor;
%             case 'a3'
%                 epsi.a3 = (epsi.a3-acc_offset)/acc_factor;
%         end
    end
end

%% give a name to field

epsi_fields = fieldnames(epsi);
for i  = 1:numel(epsi_fields)
    epsi.(epsi_fields{i}) = reshape(epsi.(epsi_fields{i}),[],1);
end

epsi.epsitime=nan.*epsi.a3;
dt=diff(efe.time);
mdt=dt./efe.data.nblocks;
try
    dt=[dt; dt(end)];
    mdt=[mdt; mdt(end)];
    for t=1:length(efe.time)
        epsi.epsitime(1+(t-1)*efe.data.nblocks: ...
            (t-1)*efe.data.nblocks+efe.data.nblocks)=linspace(efe.time(t),efe.time(t)+dt(t)-mdt(t),efe.data.nblocks);
    end
catch
    warning('no dt yet')
end


%% SBE 49
[ind_sbe_start, ind_sbe_end,ind_sbe_tokens] = regexp(str,'\$SBE([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end','tokenExtents');
