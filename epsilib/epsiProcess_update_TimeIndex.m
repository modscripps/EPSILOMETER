function [] = epsiProcess_update_TimeIndex(dirname,filename,data)

datatime = data.time_s;
if isfield(data,'dnum')
datadnum = data.dnum;
end

if exist(fullfile(dirname, 'TimeIndex.mat'),'file')
    load(fullfile(dirname, 'TimeIndex.mat'));
    ind = strncmp(filename,TimeIndex.filenames,length(filename));
    if sum(ind) ~= 1
        TimeIndex.filenames = [TimeIndex.filenames; {filename}];
        TimeIndex.timeStart = cat(1,TimeIndex.timeStart,datatime(1));
        TimeIndex.timeEnd = cat(1,TimeIndex.timeEnd,datatime(end));
        if exist('datadnum')
            TimeIndex.dnumStart = cat(1,TimeIndex.dnumStart,datadnum(1));
            TimeIndex.dnumEnd = cat(1,TimeIndex.dnumEnd,datadnum(end));
        end
    else
        TimeIndex.timeStart(ind) = datatime(1);
        TimeIndex.timeEnd(ind) = datatime(end);
        if exist('datadnum')
        TimeIndex.dnumStart(ind) = datadnum(1);
        TimeIndex.dnumEnd(ind) = datadnum(end);
        end
    end
else
    TimeIndex.filenames = {filename};
    TimeIndex.timeStart = datatime(1);
    TimeIndex.timeEnd = datatime(end);
    if exist('datadnum')
        TimeIndex.dnumStart = datadnum(1);
        TimeIndex.dnumEnd = datadnum(end);
    end
end
save(fullfile(dirname, 'TimeIndex.mat'),'TimeIndex');
end