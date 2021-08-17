% function process_streaming()

% when you are dome checking type:
% ctrl-C 
% Important you need to close the serial port
% delete(instrfind)



% open serial port 
listserial=seriallist;
for l=1:numel(listserial)
    fprintf('%i=%s\n',l,listserial{l})
end
id_port=input('Choose your port\nPort:');

 ser=serial(listserial{id_port},'BaudRate',115200);
%ser=serial(listserial{6},'BaudRate',115200);
fopen(ser);

%% plot stuff
addpath(genpath(fullfile('..','toolbox')));

% define some important and fixed variable
epsi.bytes_per_channel=3; % ADC is 3 bytes
epsi.nbsamples=160;% number of epsi blocks is 160 ~ 0.5 seconds
epsi.nbblock_diag=10;% 10*.5sec blocks = 5 seconds
epsi.name_length=5; % 5 bytes EPSI
epsi.finishblock=2; % 2 bytes \r\n
namechannels={'t1','t2', ...
    's1','s2', ...
    'c', ...
    'a1','a2','a3'};
countconversion={'Unipolar','Unipolar', ...
    'Unipolar','Unipolar', ...
    'Unipolar', ...
    'Unipolar','Unipolar','Unipolar'};

% accelerometer Voltage into Accelereation units (in g).
full_range = 2.5;
bit_counts = 24;
gain = 1;
acc_offset = 1.65;
acc_factor = 0.66;



length_diag=(epsi.nbblock_diag)*epsi.nbsamples; % 1760 sample
timeaxis=linspace(0,length_diag/325,length_diag);

% spectrum stuff
% sample rate channels
nb_segment=5;
FS        = 325;
tscan=5.5;
% number of samples per scan (1s) in channels
df        = 1/tscan;
f=(df:df:FS/2)'; % frequency vector for spectra
data=nan(length(namechannels),nb_segment,length_diag);
Fn    = .5*FS;  % Nyquist frequency
FR    = 2.5;    % Full range in Volts
def_noise=@(x)((FR/2^x)^2 /Fn);
Accelnoise=45e-6^2+0*f;

% now we want to open epsi_data.bin read and process the las 5 seconds
% we will plot time series and spectra with 5 seconds length.
count=0;
laps=0;
begin=0;
str=[];
while begin==0
    bytes = fscanf(ser,'%c',512);
    str=[str bytes];
    ind_madre1 = strfind(str,'$MADRE');
    if numel(ind_madre1)==2
        begin=1;
        str1=str(ind_madre1(end):end);
        blocksize=diff(ind_madre1)-1;
        disp(str1);
    end
end
str = char(ones(1,10*blocksize) * '0');
tempo_str=char(ones(1,2*blocksize) * '0');
nbbytes=numel(str1);
tempo_str(1:numel(str1))=str1;
counttoto=0;
while 1
    %     fid=fopen('../STREAMING/epsi_data.bin','r');
    %     str = fread(fid,'*char')';
    counttoto1=0;
    while nbbytes<=blocksize+1
        bytes = fscanf(ser,'%c',512);
        tempo_str(nbbytes+1:nbbytes+numel(bytes))=bytes;
        nbbytes=nbbytes+numel(bytes);
        counttoto1=counttoto1+1;
    end
    
    ind_madre1 = strfind(tempo_str,'$MADRE');
    str(end-blocksize+1:end)=tempo_str(1:ind_madre1(end)-2);
    str=circshift(str,-blocksize);
    %str(end-blocksize:end)
    ind_madre = strfind(str,'$MADRE');
    ind_epsi = strfind(str,'$EPSI');
    tempo_str(1:blocksize-1)=tempo_str(ind_madre1(end):end);
    nbbytes=nbbytes-blocksize-1;
    disp(counttoto)
    %     fclose(fid);
    count=mod(count,nb_segment)+1;
    
    
end




