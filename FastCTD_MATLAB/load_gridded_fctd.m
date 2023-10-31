function FCTD = load_gridded_fctd(data_dir,start_time,end_time)

start_time_1 = datevec(start_time);
start_time_1(5:6) = 0;

end_time_1 = datevec(end_time);
if(end_time_1(5) ~= 0 || end_time_1(6) ~= 0)
    end_time_1(4) = end_time_1(4)+1;
    end_time_1(5:6) = 0;
end

time_line_59_mins_space = datenum(start_time_1):59/60/24:datenum(end_time_1);
time_line_59_mins_space_dv = datevec(time_line_59_mins_space);
time_line_59_mins_space_dv(:,5:6) = 0;
time_line_1_hr_space = datenum(time_line_59_mins_space_dv);
time_line_1_hr_space = unique(time_line_1_hr_space);

FCTD = [];
for i_hour = 1:numel(time_line_1_hr_space)
    new_fname = datestr(time_line_1_hr_space(i_hour),'yyyymmdd_HHMMSS');
    if(~exist(fullfile(data_dir,['FCTD_' new_fname '.mat']),'file'))
        continue;
    end
    fctd = load(fullfile(data_dir,['FCTD_' new_fname '.mat']));
    fctd = fctd.FCTD;
    
    if(isempty(fctd))
        continue;
    end
    
    if(isempty(FCTD))
        FCTD = fctd;
    else
        FCTD = fctd_merge(FCTD,fctd);
    end
end

fctd_fieldnames = fieldnames(FCTD);
ind = FCTD.time < end_time & FCTD.time >= start_time;
for i = 1:numel(fctd_fieldnames)
    if(~isempty(strfind(fctd_fieldnames{i},'depth')))
        continue;
    elseif (~isempty(strfind(fctd_fieldnames{i},'tGrid')))
        continue;
    end
    FCTD.(fctd_fieldnames{i}) = FCTD.(fctd_fieldnames{i})(:,ind);
end

fctd_fieldnames = fieldnames(FCTD.tGrid);
ind = FCTD.tGrid.time < end_time & FCTD.tGrid.time >= start_time;
for i = 1:numel(fctd_fieldnames)
    if(~isempty(strfind(fctd_fieldnames{i},'depth')))
        continue;
    end
    FCTD.tGrid.(fctd_fieldnames{i}) = FCTD.tGrid.(fctd_fieldnames{i})(:,ind);
end

end

function fctd = fctd_merge(fctd1,fctd2)
fctd = fctd1;
fctd_fieldnames = fieldnames(fctd);
for i = 1:numel(fctd_fieldnames)
    if(~isempty(strfind(fctd_fieldnames{i},'depth')))
        continue;
    elseif(isstruct(fctd.(fctd_fieldnames{i})))
        fctd.(fctd_fieldnames{i}) = fctd_merge(fctd1.(fctd_fieldnames{i}),fctd2.(fctd_fieldnames{i}));
    else
        fctd.(fctd_fieldnames{i}) = [fctd1.(fctd_fieldnames{i}), fctd2.(fctd_fieldnames{i})];
    end
end
end