function [epsi,ctd]=mod_som_read_epsi_files(filename,Meta_Data)

% check if it is a single file or a directory and a set of files
if ischar(filename) % dir or file
    switch exist(filename,'file')
        case 2 % if it is a file
%             fid = fopen(filename,'r');
            [epsi,ctd] = mod_som_read_epsi_raw(filename,Meta_Data);
%             fclose(fid);
        case 7 % if it is a directory
            my_epsi_file = [];
            my_epsi_file = [my_epsi_file; dir(fullfile(filename,'*data*'))];
            if isempty(my_epsi_file)
                epsi = [];
                return
            else
                % prepare to read all files
                epsi = cell(size(my_epsi_file));
                ctd  = cell(size(my_epsi_file));
                % read the files in the directory
                for i = 1:length(my_epsi_file)
                    disp(['reading ' my_epsi_file(i).name]);
                    [epsi{i},ctd{i}] = mod_som_read_epsi_raw(fullfile(filename,my_epsi_file(i).name),Meta_Data);
                end
                % combine all files into one MET structure
                if sum(cellfun(@length,epsi))>0
                    epsi = mod_combine_epsi(epsi{:});
                end
                ctd1=ctd;
                if sum(cellfun(@length,ctd))>0
                    ctd = mod_combine_ctd(ctd{:});
                end
            end
        otherwise
            error('MATLAB:mod_read_epsi_raw:wrongFILENAME','Invalid file specified.');
    end
elseif iscellstr(filename) % cell of files
    % prepare to read all files
    epsi = cell(size(filename));
    ctd  = cell(size(filename));
    % read all files
    for i = 1:length(filename)
        disp(['reading ' filename{i}]);
        [epsi{i},ctd{i}] = mod_som_read_epsi_raw(filename{i},Meta_Data);
    end
    % combine all files into one epsi structure
    epsi = mod_combine_epsi(epsi{:});
    ctd  = mod_combine_ctd(ctd{:});
else
    
    if (filename<1)
        error('MATLAB:mod_read_epsi_raw:wrongFID','FID is invalid');
    end
    
    [epsi,ctd] = mod_som_read_epsi_raw(filename,Meta_Data);
    
    return
end

save(fullfile(Meta_Data.Epsipath,['epsi_' Meta_Data.deployment '.mat']),'epsi','-v7.3');
save(fullfile(Meta_Data.CTDpath, ['ctd_' Meta_Data.deployment '.mat']),'ctd','-v7.3');

end


function [epsi,ctd]=mod_som_read_epsi_raw(filename,Meta_Data)


fid = fopen(filename);
fseek(fid,0,1);
frewind(fid);
str = fread(fid,'*char')';
fclose(fid);

%%

% find pattern
% data format:
% $$EFE,status,timeEFEDATA*chksum\r\n
% token indices starts at ends at the square bracket
% $$EFE[       1         ] [ 2  ]\r\n
% [ind_efe_start, ind_efe_end,ind_efe_tokens] = regexp(str,'\$\$EFE([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end','tokenExtents');

% TODO make a case where there is no EFE or SBE
[ind_efe_start, ind_efe_end,ind_efe_tokens] = regexp(str,'\$EFE([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end','tokenExtents');

if numel(unique(ind_efe_end-ind_efe_start))>1
warning("l88 read EFE: Some blocks do NOT have right size.")
else
    disp("All blocks have the right size" )
end

efe.data.n_recs = numel(ind_efe_start);
efe.header_txt = '$EFE';
efe.time_offset = length(efe.header_txt)+4;
efe.timestamps_length = 16;
efe.data_offset = efe.time_offset+efe.timestamps_length;
efe.data.nchannels = Meta_Data.PROCESS.nb_channels;
efe.data.sample_freq = Meta_Data.AFE.FS;
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
    efe.time(i) = hex2dec(str(ind_efe_tokens{i}(1,1)+efe.time_offset+(0:efe.timestamps_length-1)));
    
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
            epsi.(wh_channel)=Bipolar(FR,epsi.([wh_channel '_count']));
        case {'Unipolar','unipolar'}
            epsi.(wh_channel)=Unipolar(FR,epsi.([wh_channel '_count']));
            
    end
    
    switch Meta_Data.AFE.(wh_channel).full_range
        case 1.8
            epsi.(wh_channel) = (epsi.(wh_channel)-acc_offset)/acc_factor;
    end
end

epsi_fields = fieldnames(epsi);
for i  = 1:numel(epsi_fields)
    epsi.(epsi_fields{i}) = reshape(epsi.(epsi_fields{i}),[],1);
end

epsi.epsitime=nan.*epsi.(wh_channel);
dt=diff(efe.time);
mdt=dt./efe.data.nblocks;

dt=[dt; dt(end)];
mdt=[mdt; mdt(end)];
for t=1:length(efe.time)
    epsi.epsitime(1+(t-1)*efe.data.nblocks: ...
                       (t-1)*efe.data.nblocks+efe.data.nblocks)=linspace(efe.time(t),efe.time(t)+dt(t)-mdt(t),efe.data.nblocks);
end
epsi.epsitime=epsi.epsitime./1000;

%% SBE 49
[ind_sbe_start, ~, ind_sbe_tokens] = regexp(str,'\$SBE([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n','start','end','tokenExtents');

if isempty(ind_sbe_start)
    disp('no ctd data')
    ctd=[];
else
    
    ctd.data.n_block = numel(ind_sbe_start);
    ctd.data.n_recs  = numel(ind_sbe_start)*Meta_Data.CTD.sample_per_record;
    ctd.header_txt   = '$SBE';
    ctd.header_chksum = uint8(0);
    ctd.status_offset = 1;
    ctd.time_offset = 6+1;
    ctd.timestamps_length=16;
    ctd.data_offset = ctd.time_offset+ctd.timestamps_length;
    % TODO fix the SBE name bug
    switch Meta_Data.CTD.name
        case{"SBE49","SBE"}
            ctd.data.format = 'eng';
            ctd.data.length=22;
            ctd.data.sample_freq = 16;
            ctd.data.P_raw      = NaN(ctd.data.n_recs,1);
            ctd.data.T_raw      = NaN(ctd.data.n_recs,1);
            ctd.data.S_raw      = NaN(ctd.data.n_recs,1);
            ctd.data.C_raw      = NaN(ctd.data.n_recs,1);
            ctd.data.PT_raw     = NaN(ctd.data.n_recs,1);
            ctd.cal=Meta_Data.CTD.cal;

        case{"SBE41"}
            ctd.data.format = 'PTS';
            ctd.data.length=26;
            ctd.data.sample_freq = 1;
    end
    
    
    ctd.data.sample_period = 1/ctd.data.sample_freq;

    block_time          = NaN(ctd.data.n_recs,1);
    ctd.data.time       = NaN(ctd.data.n_recs,1);
    ctd.data.status     = NaN(ctd.data.n_recs,1);
    
    ctd.data.P      = NaN(ctd.data.n_recs,1);
    ctd.data.T      = NaN(ctd.data.n_recs,1);
    ctd.data.S      = NaN(ctd.data.n_recs,1);
    ctd.data.C      = NaN(ctd.data.n_recs,1);
    ctd.data.checksum.data = uint8(zeros(ctd.data.n_recs,1));
    
    n_rec=0;
    for i = 1:ctd.data.n_block
        % might want to check the length of data too!!! right now I am skipping
        % that step
        tmp_sbe_block = str(2+ind_sbe_tokens{i}(1,1):ind_sbe_tokens{i}(1,2)); % +2 becasue now the header is SBE41 or SBE49 and not SBE
        block_time(i)  = hex2dec(tmp_sbe_block(1:ctd.timestamps_length))./1000;
        tmp_data_ctd=tmp_sbe_block(ctd.timestamps_length+1:end);
        
        for j=1:Meta_Data.CTD.sample_per_record
            n_rec=n_rec+1;
            ctd.data.time(n_rec)  = block_time(i) + (j-1)./ctd.data.sample_freq;
            rec_ctd=tmp_data_ctd((j-1)*ctd.data.length+(1:ctd.data.length));
            rec_ctd=rec_ctd(rec_ctd~=' ');
            
            switch ctd.data.format
                case 'PTS'
                    data_split   = strsplit(rec_ctd,',');
                    ctd.data.P(n_rec)        = str2double(data_split{1});
                    ctd.data.T(n_rec)        = str2double(data_split{2});
                    ctd.data.S(n_rec)        = str2double(data_split{3});
                    ctd.data.C(n_rec)        = NaN;
                case 'eng'
                    raw_sample   = rec_ctd;
                    ctd.data.T_raw(n_rec) = hex2dec(raw_sample(:,1:6));
                    ctd.data.C_raw(n_rec) = hex2dec(raw_sample(:,(1:6)+6));
                    ctd.data.P_raw(n_rec) = hex2dec(raw_sample(:,(1:6)+12));
                    ctd.data.PT_raw(n_rec) = hex2dec(raw_sample(:,(1:4)+18));
            end
        end
    end
    switch ctd.data.format
        case 'eng'
            ctd = sbe49_ascii_get_temperature(ctd);
            ctd = sbe49_ascii_get_pressure(ctd);
            ctd = sbe49_ascii_get_conductivity(ctd);
            ctd.data.S=sw_salt(ctd.data.C*10./sw_c3515,ctd.data.T,ctd.data.P);
            ctd.data.sig=sw_pden(ctd.data.S,ctd.data.T,ctd.data.P,0);
    
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
start_time_file=cellfun(@(x) x.data.time(1), varargin);
[~,I]=sort(start_time_file);
varargin=varargin(I);

% grab the data field name
tmp_ctd_fields0 = fieldnames(varargin{1}.data);
nfield=0;
for j = 1:numel(tmp_ctd_fields0)
    if numel(varargin{1}.data.(tmp_ctd_fields0{j}))==varargin{1}.data.n_recs
        nfield=nfield+1;
        ctd_fields{nfield} = tmp_ctd_fields0{j};
    end
end

% for i = 2:nargin
%     tmp_fields = fieldnames(varargin{i}.data);
%     for j = 1:numel(tmp_fields)
%         if ~ismember(tmp_fields{j},ctd_fields)
%             ctd_fields{end+1} = tmp_fields{j};
%         end
%     end
% end

for i=1:length(ctd_fields)
        evalstr = strcat('ctd.', ctd_fields{i}, '= [');
        for j=1:new_nargin
            if ~isfield(varargin{j}.data,(ctd_fields{i}))
                varargin{j}.data.(ctd_fields{i}) = NaN(size(varargin{j}.data.time));
            end
            evalstr = strcat(evalstr, 'varargin{', num2str(j), '}.data.', ctd_fields{i}, ';');
        end
        if ~isfield(varargin{nargin}.data,(ctd_fields{i}))
            varargin{nargin}.data.(ctd_fields{i}) = NaN(size(varargin{nargin}.data.time));
        end
        evalstr = strcat(evalstr, 'varargin{', num2str(nargin), '}.data.', ctd_fields{i}, '];');
        eval(evalstr);
end

end

%%
%  reads and apply calibration to the temperature data
function ctd = sbe49_ascii_get_temperature(ctd)

a0 = ctd.cal.ta0;
a1 = ctd.cal.ta1;
a2 = ctd.cal.ta2;
a3 = ctd.cal.ta3;

mv = (ctd.data.T_raw-524288)/1.6e7;
r = (mv*2.295e10 + 9.216e8)./(6.144e4-mv*5.3e5);
ctd.data.T = a0+a1*log(r)+a2*log(r).^2+a3*log(r).^3;
ctd.data.T = 1./ctd.data.T - 273.15;
return;
end

%  reads and apply calibration to the conductivity data
function ctd = sbe49_ascii_get_conductivity(ctd)
try 
g = ctd.cal.g;
h = ctd.cal.h;
i = ctd.cal.i;
j = ctd.cal.j;
tcor = ctd.cal.tcor;
pcor = ctd.cal.pcor;
catch
g = ctd.cal.cg;
h = ctd.cal.ch;
i = ctd.cal.ci;
j = ctd.cal.cj;
tcor = ctd.cal.ctcor;
pcor = ctd.cal.cpcor;
end

f = ctd.data.C_raw/256/1000;

ctd.data.C = (g+h*f.^2+i*f.^3+j*f.^4)./(1+tcor.*ctd.data.T+pcor.*ctd.data.P);

return;
end

%  reads and apply calibration to the pressure data
function ctd = sbe49_ascii_get_pressure(ctd)
% ALB 04112019 Changed ctd.cal.SBEcal. to ctd.cal.
pa0 = ctd.cal.pa0;
pa1 = ctd.cal.pa1;
pa2 = ctd.cal.pa2;
ptempa0 = ctd.cal.ptempa0;
ptempa1 = ctd.cal.ptempa1;
ptempa2 = ctd.cal.ptempa2;
ptca0 = ctd.cal.ptca0;
ptca1 = ctd.cal.ptca1;
ptca2 = ctd.cal.ptca2;
ptcb0 = ctd.cal.ptcb0;
ptcb1 = ctd.cal.ptcb1;
ptcb2 = ctd.cal.ptcb2;


y = ctd.data.PT_raw/13107;

t = ptempa0+ptempa1*y+ptempa2*y.^2;
x = ctd.data.P_raw-ptca0-ptca1*t-ptca2*t.^2;
n = x*ptcb0./(ptcb0+ptcb1*t+ptcb2*t.^2);

ctd.data.P = (pa0+pa1*n+pa2*n.^2-14.7)*0.689476;

return;
end

