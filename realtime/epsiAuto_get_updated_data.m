function [obj,tMaxNow] = epsiAuto_get_updated_data(oldData,newData,tMaxPrevious)

% TODO: Work with old data and new data

time_offset = 0; % SAN added to correct for current 12 hrs offset in
%data July 3rd.

% Load time index and get index of most recent file
obj = oldData;
load(fullfile(obj.Meta_Data.MATpath,'TimeIndex'));
[~,idxLast] = max(TimeIndex.timeEnd);

% Load most recent data file
% data = load(fullfile(obj.Meta_Data.MATpath,TimeIndex.filenames{idxLast}));
% fprintf([TimeIndex.filenames{idxLast},'\n'])

data = newData;

% SAN add offset just in case to correct for time offset
data.epsi.dnum = data.epsi.dnum+time_offset;
data.ctd.dnum = data.ctd.dnum+time_offset;
data.alt.dnum = data.alt.dnum+time_offset;
% Get last epsi dnum (TimeIndex.timeEnd is in seconds, not
% datenum, and I'd rather work in datenum)
tMaxNow.epsi = nanmax(data.epsi.dnum);
tMaxNow.ctd = nanmax(data.ctd.dnum);
tMaxNow.alt = nanmax(data.alt.dnum);

% Get indices of old data (everything up until tMaxPrevious)
idxOldEpsi = obj.epsi.dnum<=tMaxPrevious.epsi;
idxOldCtd = obj.ctd.dnum<=tMaxPrevious.ctd;
idxOldAlt = obj.alt.dnum<=tMaxPrevious.alt;
% Get length of old data (this is not the sum because there could be nans
% in the middle!)
nOldEpsi = sum(idxOldEpsi);
nOldCtd = sum(idxOldCtd);
nOldAlt = sum(idxOldAlt);

% Get indices of new data (everything after tMaxPrevious)
idxNewEpsi = data.epsi.dnum>tMaxPrevious.epsi;
idxNewCtd = data.ctd.dnum>tMaxPrevious.ctd;
idxNewAlt = data.alt.dnum>tMaxPrevious.alt;
% Get length of new data (this is not the sum because there could be nans
% in the middle!)
nNewEpsi = sum(idxNewEpsi);
nNewCtd = sum(idxNewCtd);
nNewAlt = sum(idxNewAlt);

% Get fields to put new data into structures
epsiFields = fields(obj.epsi);
ctdFields = fields(obj.ctd);
altFields = fields(obj.alt);

% fileName = ['line38' '_',...
%     num2str(hour(datetime(tMaxNow,'convertfrom','datenum')),'%02.0f'),...
%     num2str(minute(datetime(tMaxNow,'convertfrom','datenum')),'%02.0f'),...
%     num2str(second(datetime(tMaxNow,'convertfrom','datenum')),'%02.0f')];
% save(fullfile('/Volumes/Berry/blt_epsi/whats_wrong', fileName))

%% Put new epsi data into obj.epsi
for iField=1:length(epsiFields)
    % Shift old data backwards by nEpsi
    obj.epsi.(epsiFields{iField})(1:end-nNewEpsi) = obj.epsi.(epsiFields{iField})(nNewEpsi+1:end);
    
    % Where did we leave off after shifting?
    lastIdx = nOldEpsi-nNewEpsi;
    if lastIdx<0 %If lastIdx is negative, it's the first time you're calling it
        lastIdx = 0;
    end
    
    % Now add the new data
    obj.epsi.(epsiFields{iField})(lastIdx+1:lastIdx+nNewEpsi) = data.epsi.(epsiFields{iField})(idxNewEpsi);
end

%% Put new ctd data into obj.ctd
for iField=1:length(ctdFields)
    % Shift old data backwards by nCtd
    obj.ctd.(ctdFields{iField})(1:end-nNewCtd) = obj.ctd.(ctdFields{iField})(nNewCtd+1:end);
    
    % Where did we leave off after shifting?
    lastIdx = nOldCtd-nNewCtd;
    if lastIdx<0 %If lastIdx is negative, it's the first time you're calling it
        lastIdx =0;
    end
    
    % Now add the new data
    obj.ctd.(ctdFields{iField})(lastIdx+1:lastIdx+nNewCtd) = data.ctd.(ctdFields{iField})(idxNewCtd);
end

%% Put new alt data into obj.alt
for iField=1:length(altFields)
    % Shift old data backwards by nAlt
    obj.alt.(altFields{iField})(1:end-nNewAlt) = obj.alt.(altFields{iField})(nNewAlt+1:end);
       
    % Where did we leave off after shifting?
    lastIdx = nOldAlt-nNewAlt;
    if lastIdx<0 %If lastIdx is negative, it's the first time you're calling it
        lastIdx = 0;
    end
    
    % Now add the new data
    obj.alt.(altFields{iField})(lastIdx+1:lastIdx+nNewAlt) = data.alt.(altFields{iField})(idxNewAlt);
end


clear data
% 
% fileName = ['line102' '_',...
%     num2str(hour(datetime(tMaxNow,'convertfrom','datenum')),'%02.0f'),...
%     num2str(minute(datetime(tMaxNow,'convertfrom','datenum')),'%02.0f'),...
%     num2str(second(datetime(tMaxNow,'convertfrom','datenum')),'%02.0f')];
% save(fullfile('/Volumes/Berry/blt_epsi/whats_wrong', fileName))

end