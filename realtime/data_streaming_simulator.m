function [] = data_streaming_simulator(file_dir_in,file_dir_out,varargin)

% Simulate streaming epsi data by slowly copying the contents of one raw
% file directory into a new one. This way, you can test realtime scripts on
% old data sets without needing to have an epsi actively collecting data.
%
% INPUTS
%   file_dir_in - input directory where existing raw epsi files live
%   file_dir_out - output directory where you want to copy files simulating
%                  data coming in
%   suffixSearch (optional) - character string to look for in raw files
%                             (only copy the files that contain that string)
%
% Nicole Couto | February 2022
% -------------------------------------------------------------------------
% Change these parameters to adjust how quickly data "streams" in.
%   n_blocks:           number of lines that get copied before pausing
%   seconds_between_blocks: amount of time in seconds to pause before copying
%                           the next block
n_blocks = 100; %Number of blocks that will be copied before pausing
seconds_between_blocks = 15; %Number of seconds to pause between blocks
do_pause = 1; %Switch pause on/off
str_to_find = {'$EFE','$SB49','$SB41','$ALT'}; %Add more data header strings here if you want to copy them

if ~exist(file_dir_out,'dir')
    eval([ '!mkdir ' strrep(file_dir_out,' ','\ ')]);
end

% Check for suffix. If nothing is specified, use all files
if nargin<5
    suffixSearch = '*.*';
    suffixStr = '.';
else
    suffixStr = varargin{1};
    suffixSearch = ['*',suffixSearch];
end

%% List the ascii files in the directory
myASCIIfiles = dir(fullfile(file_dir_in, suffixSearch));

%% EFE offset values
tag.sync.strvalue='$';
tag.sync.offset=0;
tag.sync.length=1;

tag.header.strvalue  = 'FFFF';
tag.header.length = strlength(tag.header.strvalue);
tag.header.offset = tag.sync.offset+tag.sync.length;

tag.hextimestamp.strvalue  = "0000000000000000";
tag.hextimestamp.length = strlength(tag.hextimestamp.strvalue);
tag.hextimestamp.offset = tag.header.offset+tag.header.length;

tag.hexlengthblock.strvalue  = "00000000";
tag.hexlengthblock.length = strlength(tag.hexlengthblock.strvalue);
tag.hexlengthblock.offset = tag.hextimestamp.offset+tag.hextimestamp.length;

tag.headerchecksum.strvalue = "*FF";
tag.headerchecksum.length   = strlength(tag.headerchecksum.strvalue);
tag.headerchecksum.offset   = tag.hexlengthblock.offset+tag.hexlengthblock.length;

tag.data_offset = tag.headerchecksum.offset+tag.headerchecksum.length+1;
tag.chksum.strvalue = "FFFFF" ;
tag.chksum.length   = strlength(tag.chksum.strvalue);

efe.data.n_channels         = 7;
efe.data.sample_freq        = 320;
efe.data.sample_period      = 1/efe.data.sample_freq;
efe.data.bytes_per_channel  = 3;
efe.data.timestamp_length   = 8;
efe.data.n_elements         = efe.data.timestamp_length + ...
    efe.data.n_channels*efe.data.bytes_per_channel; % 8 bytes timestamps + 3 bytes ADC
efe.data.recs_per_block     = 80;

T_length                    = 11; %Length of T timestamp header before all blocks

%% Loop through files and copy contents into new files in a file_dir_out
for idx=1:length(myASCIIfiles)

    % Disregard the file if it starts with '.'
    if strcmp(myASCIIfiles(idx).name(1),'.')
        continue
    end
    
    myRAWfile = fullfile(file_dir_in,myASCIIfiles(idx).name);
    fid = fopen(myRAWfile,'rb');

    % Read the entire file as in mod_som_read_epsi_files.m. Find the
    % indices of every data block
    if fid>0
        fseek(fid,0,1);
        frewind(fid);
        str = fread(fid,'*char')';

        % Find the start and stop indices for every block type. Concatenate them
        % all so you can sort and write blocks in order of their start index.
        blocks.ind_start = [];
        blocks.ind_stop = [];
        s1 = 1;
        for ii=1:length(str_to_find)
            [ind_start, ind_stop]   = regexp(str,['\'  str_to_find{ii} '([\S\s]+?)\*([0-9A-Fa-f][0-9A-Fa-f])\r\n'],'start','end');
            if ~isempty(ind_start) && length(ind_start)==length(ind_stop)
                blocks.ind_start = [blocks.ind_start,ind_start];
                blocks.ind_stop = [blocks.ind_stop,ind_stop];

                % Copy the header string to every index. Then adjust position s1.
                s2 = s1+length(blocks.ind_start)-1;
                [blocks.header_str{s1:s2}] = deal(str_to_find{ii});
                s1 = s2+1;
            end
        end

        if ~isempty(ind_start)
            % Sort the indices and associated header strings by ind_start
            [~,iSort] = sort(blocks.ind_start);
            sorted_blocks.ind_start = blocks.ind_start(iSort);
            sorted_blocks.ind_stop = blocks.ind_stop(iSort);
            sorted_blocks.header_str = blocks.header_str(iSort);

            % Print every data block in order. Use fprintf for non-epsi and fwrite
            % for binary epsi? Arnaud suggests 'char' or 'uint8'. Perhaps I can use
            % fwrite for every line if char works.

            % Initialize a line counter
            line_count = 0;

            % Create the new raw file
            newRAWfile = fullfile(file_dir_out,myASCIIfiles(idx).name);
            new_fid = fopen(newRAWfile,'wb');
            % loop until end of file is reached
            for b=1:length(sorted_blocks.ind_start)

                if ~strcmp(sorted_blocks.header_str{b},'$EFE')
                    % If data is non-binary, print it to the new file with
                    % fprintf
                    fprintf(new_fid,'\n%s',str(sorted_blocks.ind_start(b)-T_length:sorted_blocks.ind_stop(b)));

                elseif strcmp(sorted_blocks.header_str{b},'$EFE')
                    % If data is binary, break it up into header and data.
                    % Then break up the data into bytes and print it byte
                    % by byte

                    % Grab the block of data starting with the header
                    efe_block_str = str(sorted_blocks.ind_start(b):sorted_blocks.ind_stop(b));

                    % Get the data after the header.
                    efe_block_data = efe_block_str(tag.data_offset:end-tag.chksum.length);

                    % First, write the header
                    fprintf(new_fid,'\n%s',...
                        str(sorted_blocks.ind_start(b)-T_length:sorted_blocks.ind_start(b)+tag.data_offset-2));
                    % Then write the timestamp and ADC data
                    fwrite(new_fid,uint32(efe_block_data),'uchar');
                    % Finally, add the checksum at the end
                    fprintf(new_fid,'%s',efe_block_str(end-4:end));

                end

                line_count = line_count+1;

                if mod(line_count,n_blocks)==0 && do_pause
                    pause(seconds_between_blocks);
                end
            end

            % close the files
            fclose(fid);
            fclose(new_fid);

        end
    end

end
