function obj = processNewProfiles(obj,varargin)
% obj = f_processNewProfiles(obj,varargin)
%
% INPUTS:
%   obj = epsi_class object
%
% OPTIONAL INPUTS:
%   'grid',P - 'grid' flag to grid some of the profile variables onto a
%               standard pressure grid 
%            - P = the pressure array to use
%
% If there are optional arguments they are 'grid' and P, the pressure array
% to grid onto
if nargin>1
    makeGrid = true;
    P = varargin{1}{2};
else
    makeGrid = false;
end

% % Pick out profile drops from CTD pressure timeseries
% obj.f_getProfileIndices;

% Load PressureTimeseries
load(fullfile(obj.Meta_Data.MATpath,'PressureTimeseries.mat'))

% Look for the current list of profiles
profList = dir(fullfile(obj.Meta_Data.L1path,'Profile*.mat'));
profNumChar = cell2mat(cellfun(@(C) C(8:10),{profList(:).name},'uniformoutput',0).');
if ~isempty(profNumChar)
    lastProfNum = str2double(profNumChar(end,:));
    % Load the last profile
    lastProf = load(fullfile(obj.Meta_Data.L1path,sprintf('Profile%03.f',lastProfNum)));
    
    % Does the last profile have all its data? Or was more collected in the
    % last batch of files?
    lastProf_maxTime = nanmax(lastProf.Profile.ctddnum);
    
    % If there is more data for the last profile in the PressureTimeseries,
    % rerun the last profile. In the next step, we'll run everything after the
    % lastProfNum, so subtract one from that value.
    if PressureTimeseries.ctddnum(PressureTimeseries.endprof(lastProfNum))>lastProf_maxTime
        lastProfNum = lastProfNum-1;
    end
elseif isempty(profNumChar)
    lastProfNum = 0;
end


% Loop through the profile indices in PressureTimeseries. Process the new
% ones
for iProf=1:length(PressureTimeseries.startprof)
    if iProf>lastProfNum
            profIdx = PressureTimeseries.startprof(iProf):PressureTimeseries.endprof(iProf);
            tMin = PressureTimeseries.ctddnum(PressureTimeseries.startprof(iProf));
            tMax = PressureTimeseries.ctddnum(profIdx(end));
            fprintf('Building Profile%03.0f of %03.0f\n',iProf,length(PressureTimeseries.startprof))
            
            Profile = obj.f_cropTimeseries(tMin,tMax);
            Profile.profNum = iProf;
            Profile = obj.f_computeTurbulence(Profile);
            
            if makeGrid
                process_gridProfiles;
            end
            
            % Save new profile
            saveName = fullfile(obj.Meta_Data.L1path,sprintf('Profile%03.0f',iProf));
            eval(['save ' saveName ' Profile']);
            clear Profile
    end
end
end %end f_processNewProfiles
