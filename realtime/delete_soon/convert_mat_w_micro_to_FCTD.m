function [] = convert_mat_w_micro_to_FCTD(dirs,Meta_Data,input_struct)
% dirs.mat
%     .fctd_mat

% Compare TimeIndex.mat to FastCTD_MATfile_TimeIndex.mat to determine which
% files to convert. Always convert the last file in
% FastCTD_MATfile_TimeIndex, because it may have been created from a .mat
% file that hasn't finished being filled. If there are no FastCTD files
% yet, process all the files in the .mat directory.
load(fullfile(dirs.mat,'TimeIndex.mat'));
if exist(fullfile(dirs.fctd_mat,'FastCTD_MATfile_TimeIndex.mat'),'file')
    load(fullfile(dirs.fctd_mat,'FastCTD_MATfile_TimeIndex.mat'));
else
    FastCTD_MATfile_TimeIndex.filenames = {''};
end

process_files = setdiff(TimeIndex.filenames,FastCTD_MATfile_TimeIndex.filenames);
if ~isempty(process_files)
    process_files = TimeIndex.filenames;
else
    process_files = [process_files;FastCTD_MATfile_TimeIndex.filenames{end}];
end

if ~isempty(process_files)
    for i=1:length(process_files)
        matData = load(fullfile(dirs.mat,process_files{i}));
        base = process_files{i};

        % Calculate turbulence in .mat files
            matData.micro = epsiProcess_calc_turbulence(Meta_Data,matData,0);
            save(fullfile(dirs.mat,base),'-struct','matData')

        % Make FCTD mat files, including epsilon and chi data
        make_FCTD_mat(matData,dirs.fctd_mat,base,input_struct.cruise_specifics);
    end
end