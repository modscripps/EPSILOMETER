function [obj,tMaxNow] = getUpdatedData(obj,tMaxPrevious)

time_offset = 0; % SAN added to correct for current 12 hrs offset in data

% Load time index and get most recent time value
load(fullfile(obj.Meta_Data.MATpath,'Epsi_MATfile_TimeIndex'));
[~,idxLast] = max(Epsi_MATfile_TimeIndex.timeEnd);

% Load most recent data file
data = load(fullfile(obj.Meta_Data.MATpath,Epsi_MATfile_TimeIndex.filenames{idxLast}));
fprintf([Epsi_MATfile_TimeIndex.filenames{idxLast},'\n'])


% SAN add offset just in case to correct for time offset
data.epsi.epsidnum = data.epsi.epsidnum+time_offset;
data.ctd.ctddnum = data.ctd.ctddnum+time_offset;
data.alt.altdnum = data.alt.altdnum+time_offset;

% Get last epsi dnum (Epsi_MATfile_TimeIndex.timeEnd is in seconds, not
% datenum, and I'd rather work in datenum)
tMaxNow = nanmax(data.epsi.epsidnum);

% Get indices of new data (everything between tMaxPrevious and tMaxNow)
idxDataEpsi = data.epsi.epsidnum>tMaxPrevious;
idxDataCtd = data.ctd.ctddnum>tMaxPrevious;
idxDataAlt = data.alt.altdnum>tMaxPrevious;

% Get length of new data
nEpsi = sum(idxDataEpsi);
nCtd = sum(idxDataCtd);
nAlt = sum(idxDataAlt); 

% Get fields to put new data into structures
epsiFields = fields(obj.epsi);
ctdFields = fields(obj.ctd);
altFields = fields(obj.alt);

%% Put new epsi data into obj.epsi
for iField=1:length(epsiFields)
    % Shift old data backwards by nEpsi
   obj.epsi.(epsiFields{iField})(1:nEpsi) = obj.epsi.(epsiFields{iField})(nEpsi+1:2*nEpsi);
   
%    % Make everything after it nan. This is important if the new data chunk is
%    % smaller than the old data chunk. After shift, you would have left some
%    % lingering old data at the end. Get rid of it!
%    obj.epsi.(epsiFields{iField})((2*nEpsi)+nEpsi:end) = nan;
   
   % Now add the new data
   obj.epsi.(epsiFields{iField})((2*nEpsi)+1:(2*nEpsi)+nEpsi) = data.epsi.(epsiFields{iField})(idxDataEpsi);
end

% There might have been some bad indexing, so sort the date by datenum. 
[~,idxOrder] = sort(obj.epsi.epsidnum);
for iField=1:length(epsiFields)
   obj.epsi.(epsiFields{iField}) = obj.epsi.(epsiFields{iField})(idxOrder);
end

%% Put new ctd data into obj.ctd
for iField=1:length(ctdFields)
    % Shift old data backwards by nCtd
   obj.ctd.(ctdFields{iField})(1:nCtd) = obj.ctd.(ctdFields{iField})(nCtd+1:2*nCtd);
   
%    % Make everything after it nan. This is important if the new data chunk is
%    % smaller than the old data chunk. After shift, you would have left some
%    % lingering old data at the end. Get rid of it!
%    obj.ctd.(ctdFields{iField})((2*nCtd)+nCtd :end) = nan;
    
   % Now add the new data
   obj.ctd.(ctdFields{iField})((2*nCtd)+1:(2*nCtd)+nCtd) = data.ctd.(ctdFields{iField})(idxDataCtd);
end
% There might have been some bad indexing, so sort the date by datenum. 
[~,idxOrder] = sort(obj.ctd.ctddnum);
for iField=1:length(ctdFields)
   obj.ctd.(ctdFields{iField}) = obj.ctd.(ctdFields{iField})(idxOrder);
end

%% Put new alt data into obj.alt
for iField=1:length(altFields)
    % Shift old data backwards by nAlt
   obj.alt.(altFields{iField})(1:nAlt) = obj.alt.(altFields{iField})(nAlt+1:2*nAlt);
   
%    % Make everything after it nan. This is important if the new data chunk is
%    % smaller than the old data chunk. After shift, you would have left some
%    % lingering old data at the end. Get rid of it!
%    obj.alt.(altFields{iField})((2*nAlt)+nAlt:end) = nan;
   
   % Now add the new data
   obj.alt.(altFields{iField})((2*nAlt)+1:(2*nAlt)+nAlt) = data.alt.(altFields{iField})(idxDataAlt);
end
% There might have been some bad indexing, so sort the date by datenum. 
[~,idxOrder] = sort(obj.alt.altdnum);
for iField=1:length(altFields)
   obj.alt.(altFields{iField}) = obj.alt.(altFields{iField})(idxOrder);
end

clear data

end