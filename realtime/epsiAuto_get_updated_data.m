function [obj,tMaxNow] = epsiAuto_get_updated_data(oldData,newData,tMaxPrevious)

if isempty(newData.epsi)
    disp('There is no new epsi data.')
    obj = oldData;
    tMaxNow = tMaxPrevious;
else
    % Determine whether tMaxPrevious is a datenum or seconds since power on.
    % Depending on what you have, use either 'time_s' or 'dnum' as the
    % timestamp.
    if tMaxPrevious.epsi>7e5
        timestamp = 'dnum';
    else
        timestamp = 'time_s';
    end
    
    time_offset = 0; % SAN added to correct for current 12 hrs offset in
    %data July 3rd.
    
    obj = oldData;
    data = newData;
    
    % Load most recent data file
    % data = load(fullfile(obj.Meta_Data.paths.mat_data,TimeIndex.filenames{idxLast}));
    % fprintf([TimeIndex.filenames{idxLast},'\n'])
    
    
    % TODO 08/18/2021
    % We are dealing with each peripheral individually.
    % In the futur we should loop over the fileds inside the data struct
    % The new variable names (e.g., idxOldEpsi) should changes too to reflect
    % that change (e.g., idxOld.epsi).
    periphNames = {'epsi','ctd','alt'};
    for p=1:length(periphNames)
        periph = periphNames{p};
        if ~isempty(data.(periph))
            % do all the things
        end
    end
    
    % SAN add offset just in case to correct for time offset
    if ~isempty(data.epsi)
        data.epsi.(timestamp) = data.epsi.(timestamp)+time_offset;
        % Get last epsi timestamp
        tMaxNow.epsi = nanmax(data.epsi.(timestamp));
        % Get indices of old data (everything up until tMaxPrevious)
        idxOldEpsi = obj.epsi.(timestamp)<=tMaxPrevious.epsi;
        % Get length of old data (this is not the sum because there could be nans
        % in the middle!)
        nOldEpsi = sum(idxOldEpsi);
        % Get indices of new data (everything after tMaxPrevious)
        idxNewEpsi = data.epsi.(timestamp)>tMaxPrevious.epsi;
        % Get length of new data (this is not the sum because there could be nans
        % in the middle!)
        nNewEpsi = sum(idxNewEpsi);
        % Get fields to put new data into structures
        epsiFields = fields(obj.epsi);
    end
    
    if ~isempty(data.ctd)
        data.ctd.(timestamp) = data.ctd.(timestamp)+time_offset;
        % Get last ctd timestamp
        tMaxNow.ctd = nanmax(data.ctd.(timestamp));
        % Get indices of old data (everything up until tMaxPrevious)
        idxOldCtd = obj.ctd.(timestamp)<=tMaxPrevious.ctd;
        % Get length of old data (this is not the sum because there could be nans
        % in the middle!)
        nOldCtd = sum(idxOldCtd);
        % Get indices of new data (everything after tMaxPrevious)
        idxNewCtd = data.ctd.(timestamp)>tMaxPrevious.ctd;
        % Get length of new data (this is not the sum because there could be nans
        % in the middle!)
        nNewCtd = sum(idxNewCtd);
        % Get fields to put new data into structures
        ctdFields = fields(obj.ctd);
    end
    
    if ~isempty(data.alt)
        data.alt.(timestamp) = data.alt.(timestamp)+time_offset;
        % Get last alt timestamp
        tMaxNow.alt = nanmax(data.alt.(timestamp));
        % Get indices of old data (everything up until tMaxPrevious)
        idxOldAlt = obj.alt.(timestamp)<=tMaxPrevious.alt;
        % Get length of old data (this is not the sum because there could be nans
        % in the middle!)
        nOldAlt = sum(idxOldAlt);
        % Get indices of new data (everything after tMaxPrevious)
        idxNewAlt = data.alt.(timestamp)>tMaxPrevious.alt;
        % Get length of new data (this is not the sum because there could be nans
        % in the middle!)
        nNewAlt = sum(idxNewAlt);
        % Get fields to put new data into structures
        altFields = fields(obj.alt);
    end
    
    
    %% Put new epsi data into obj.epsi
    if ~isempty(data.epsi)
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
    end
    %% Put new ctd data into obj.ctd
    if ~isempty(data.ctd)
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
    end
    
    %% Put new alt data into obj.alt
    if ~isempty(data.alt)
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
    end
    
    clear data
    
end %if isempty(newData.epsi)
end