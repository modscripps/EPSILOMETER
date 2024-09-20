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

    obj = oldData;

    % Get sig
    if isfield(newData.ctd,'sgth')
        newData.ctd.sig=newData.ctd.sgth;
    end
    data = newData;

    % List data fields to add obj
    periphNames = {'epsi','ctd','alt','vnav','gps','ttv'};

    %% Get indices and lengths of old and new data for all peripherals
    for p=1:length(periphNames)
        periph = periphNames{p};
        if ~isempty(data.(periph))
            % Get length of array (this should always stay the same)
            nArray.(periph) = numel(obj.(periph).(timestamp));
            % Get last timestamp
            tMaxNow.(periph) = nanmax(data.(periph).(timestamp));
            % Get indices of old data (everything up until tMaxPrevious)
            idxOld.(periph) = obj.(periph).(timestamp)<=tMaxPrevious.(periph);
            % Get length of old data
            nOld.(periph) = sum(idxOld.(periph));
            % Get indices of new data (everything after tMaxPrevious)
            idxNew.(periph) = data.(periph).(timestamp)>tMaxPrevious.(periph);
            % Get length of new data
            nNew.(periph) = sum(idxNew.(periph));
            % Get fields to put new data into old structures
            field_list.(periph) = fields(obj.(periph));
        end
    end

    %% ALB if c1 f1 volt exist = FCTD
    if find(cellfun(@(x) contains(x,'c1_volt'),fields(data.epsi)))
        %ALB TODO figure out a way to keep c1_volt and f1_volt 
        % This does not seem to bother the blue matlab 
        % (after minor changes in make_FCTD_mat .line 107)
        data.epsi.s1_volt=data.epsi.f1_volt;
        data.epsi.s2_volt=data.epsi.c1_volt;
    end

    %% Put new data from all peripherals into obj
    for p=1:length(periphNames)
        periph = periphNames{p};
        if ~isempty(data.(periph))
            for iField=1:length(field_list.(periph))

                % Append new data to old data
                obj.(periph).(field_list.(periph){iField})(1:nOld.(periph)+nNew.(periph),:) ...
                            = [obj.(periph).(field_list.(periph){iField})(idxOld.(periph),:);...
                               data.(periph).(field_list.(periph){iField})(idxNew.(periph),:)];

                % Only keep the last nArray indices
                obj.(periph).(field_list.(periph){iField}) = ...
                    obj.(periph).(field_list.(periph){iField})(end-nArray.(periph)+1:end,:);

            end %end loop through data fields
        end %end if that periph exists
    end %end loop through periphs

end %if isempty(newData.epsi)
