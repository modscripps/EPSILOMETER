function Meta_Data= read_MetaProcess(Meta_Data,filename)
% Set temp probe numbers
fid=fopen(filename);
str_meta_data_process=textscan(fid,'%s','Delimiter','\n');
str_meta_data_process=[str_meta_data_process{:}];
fclose(fid);
for f=1:length(str_meta_data_process)
    if~isempty(str_meta_data_process{f}) && isempty(strfind(str_meta_data_process{f},'%'))
        eval(['Meta_Data.PROCESS.' str_meta_data_process{f}])
    end
end

Meta_Data.PROCESS.filename = filename;

if isfield(Meta_Data.paths,'data')
    save(fullfile(Meta_Data.paths.data,'Meta_Data.mat'),'Meta_Data')
end