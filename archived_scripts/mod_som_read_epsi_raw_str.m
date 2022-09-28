function epsi=mod_som_read_epsi_raw_str(str,Meta_Data)
%%

% find pattern
% data format:
% $$EFE,status,timeEFEDATA*chksum\r\n
% token indices starts at ends at the square bracket
% $$EFE[       1         ] [ 2  ]\r\n
[ind_efe_start, ind_efe_end,ind_efe_tokens] = regexp(str,'\$EFErev3([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end','tokenExtents');

efe.header_txt = '$EFErev3';
efe.data.n_recs = numel(ind_efe_start);
efe.status_offset = 1;
efe.time_offset = 4;
efe.data_offset = efe.time_offset+16;
efe.data.nchannels = 7;
efe.data.sample_freq = 320;
efe.data.sample_period = 1/efe.data.sample_freq;
efe.data.bytes_per_channel = 3;
efe.data.nblocks = 160;
efe.data.length = efe.data.nblocks*efe.data.nchannels*efe.data.bytes_per_channel;
efe.data.raw_bytes = int32(zeros(efe.data.n_recs,efe.data.length));
efe.n_recs = numel(ind_efe_tokens);
efe.time = NaN(efe.data.n_recs,1);
efe.status = NaN(efe.data.n_recs,1);
efe.checksum.data = uint8(zeros(efe.data.n_recs,1));
efe.header_chksum = uint8(0);

% calculate 
for i = 1:numel(efe.header_txt)
    efe.header_chksum = bitxor(efe.header_chksum,uint8(efe.header_txt(i)));
end
efe.checksum.calculate = uint8(ones(efe.data.n_recs,1))*efe.header_chksum;
% efe.checksum.calculate = uint8(zeros(efe.data.n_recs,1));
% efe.checksum.data = efe.checksum.calculate;

for i = 1:efe.data.n_recs
    % might want to check the length of data too!!! right now I am skipping
    % that step
    tmp_data_u8 = uint8(str(ind_efe_tokens{i}(1,1):ind_efe_tokens{i}(1,2)));
%     tmp_data_u8 = uint8(str(ind_efe_start(i)+1:(ind_efe_end(i))-5));
    try
        efe.data.raw_bytes(i,:) = int32(str((ind_efe_tokens{i}(1,1)+efe.data_offset):ind_efe_tokens{i}(1,2)));
    catch
        efe.data.raw_bytes(i,:) = nan;
        fprintf('jump in the stream, block %i\n',i)
    end
    efe.checksum.data(i,:) = hex2dec(str(ind_efe_tokens{i}(2,1):ind_efe_tokens{i}(2,2)));
    for j = 1:numel(tmp_data_u8) %% check sum is not working yet
        efe.checksum.calculate = bitxor(efe.checksum.calculate,tmp_data_u8(j));
    end
    efe.time(i) = hex2dec(str(ind_efe_tokens{i}(1,1)+efe.time_offset+(0:15)));
    
end


% epsi.raw1 = epsi.raw(:,1:epsi.bytes_per_channel:end)*256^2+ ...
%             epsi.raw(:,2:epsi.bytes_per_channel:end)*256+ ...
%             epsi.raw(:,3:epsi.bytes_per_channel:end);

efe.data.raw_24bits = efe.data.raw_bytes(:,1:3:end)*256^2+ ...
                      efe.data.raw_bytes(:,2:3:end)*256+ ...
                      efe.data.raw_bytes(:,3:3:end);

efe.data.raw = uint32(NaN(efe.data.n_recs,efe.data.nblocks,efe.data.nchannels));

for i = 1:efe.data.nchannels
    efe.data.raw(:,:,i) = efe.data.raw_24bits(:,i:efe.data.nchannels:end);
end

[N,M,C]=size(efe.data.raw);
efe.data.raw=permute(efe.data.raw,[2 1 3]);
efe.data.raw=reshape(efe.data.raw,[N*M C]);

%% convert count 2 volt

bit_counts = 24;
gain = 1;
acc_offset = 1.8/2;
acc_factor = 0.5;
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
