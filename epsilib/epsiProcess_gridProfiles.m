function obj = epsiProcess_gridProfiles(obj,z)

% NC 9/9/21 To do - This could be faster. Right now it always interpolates
% every profile and concatenates them. Could check to see if current
% profile already exists in grid

% Try loading griddedProfiles. If it doesn't exist, we'll
% make it
if exist(fullfile(obj.Meta_Data.paths.profiles,'griddedProfiles.mat'),'file')==2
    grid_exists = 1;
    G = load(fullfile(obj.Meta_Data.paths.profiles,'griddedProfiles'));
    GRID = G.GRID;
else
    grid_exists = 0;
end
% Actually, always re-grid from the beginning since you might be
% changing the depth array
grid_exists = 0;
clear GRID

fileList = dir(fullfile(obj.Meta_Data.paths.profiles,'Profile*.mat'));
for iFile=1:length(fileList)

    disp(sprintf('Gridding %03.0f of %03.f',iFile,length(fileList)))

    load(fullfile(obj.Meta_Data.paths.profiles,fileList(iFile).name))

    if isfield(Profile,'pr') && ~isempty(Profile.pr)%Only add to

        % Interpolate this profile to standard pressure grid
        gridNew = epsiProcess_interpolate_Profile_to_P(Profile,z);

        % Add gridded profile to full GRID. If full grid doesn't
        % exist yet, create it.
        if grid_exists
            varList = fields(GRID);
            for iVar=1:length(varList)
                if ~strcmp(varList{iVar},'mission') && ~strcmp(varList{iVar},'vehicle_name') && ~strcmp(varList{iVar},'deployment') && ~strcmp(varList{iVar},'filenames')
                    GRID.(varList{iVar})(:,end+1) = gridNew.(varList{iVar});
                end

            end
        else
            varList = fields(gridNew);
            for iVar=1:length(varList)
                GRID.(varList{iVar})(:,1) = gridNew.(varList{iVar});
            end
            grid_exists = 1;
        end

        % Keep only unique profiles - if you're running this in
        % realtime you could have part of a profiles followed by the
        % rest of it.
        % ... but keep the last one if there is more than one.
        % TO DO: make this step less janky
        profNumX = fliplr(1:length(GRID.profNum));
        profNumY = fliplr(GRID.profNum);
        [~,u] = unique(profNumY);
        for iVar=1:length(varList)
            if ~strcmp(varList{iVar},'mission') && ~strcmp(varList{iVar},'vehicle_name') ...
                    && ~strcmp(varList{iVar},'deployment') && ~strcmp(varList{iVar},'pr') ...
                    && ~strcmp(varList{iVar},'z') && ~strcmp(varList{iVar},'filenames')
                GRID.(varList{iVar}) = GRID.(varList{iVar})(:,profNumX(u));
            end
        end

        close

        GRID.pr = GRID.pr(:,1);
        GRID.z = GRID.z(:,1);
        GRID.mission = GRID.mission(:).';
        GRID.vehicle_name = GRID.vehicle_name(:).';
        GRID.deployment = GRID.deployment(:).';
        GRID.filenames{iFile} = Profile.filenames;


%{
 %% Add Meta_Data
GRID.Meta_Data = Meta_Data;

%% Add varInfo
GRID.varInfo.pr = {'CTD pressure','db'};
GRID.varInfo.dnum = {'datenum','Matlab datenum'};
GRID.varInfo.w = {'fall speed','db s^{-1}'};
GRID.varInfo.t = {'temperature','C'};
GRID.varInfo.s = {'salinity','psu'};
GRID.varInfo.kvis = {'kinematic viscosity',''};
GRID.varInfo.epsilon = {'turbulent kinetic energy dissipation rate calculated from Ps_shear_k', ''};
GRID.varInfo.epsilon_co = {'turbulent kinetic energy dissipation rate calculated from Ps_shear_co_k', ''};
GRID.varInfo.chi = {'temperature gradient dissipation rate','°C^2 s^{-1}'};
GRID.varInfo.sh_fc = {'shear cutoff frequency, 1=uncorrected, 2=coherence-corrected', 'Hz'};
GRID.varInfo.tg_fc = {'temperature gradient cutoff frequency, 1=uncorrected, 2=coherence-corrected', 'Hz'};
GRID.varInfo.flag_tg_fc = {'temperature gradient cut off frequency is very high','0/1'};
GRID.varInfo.Coh_10_45Hz_avg = {'average coherences between 10 and 45 Hz',''};
GRID.varInfo.Pa_g_10_45Hz_avg = {'average acceleration power spectral densities between 10 and 45 Hz',''};
GRID.varInfo.Ps_volt_10_45Hz_avg = {'average shear power spectral densities between 10 and 45 Hz',''};
GRID.varInfo.Pt_volt_10_45Hz_avg = {'average temperature power spectral densities between 10 and 45 Hz',''};
GRID.varInfo.iScan = {'scan indices',''};
 
%}

%% Save data   
    saveName = fullfile(obj.Meta_Data.paths.profiles,'griddedProfiles.mat');
    save(saveName,'GRID');

    end %End if there is Profile.pr field

end %End loop through profiles