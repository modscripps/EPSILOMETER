% Add Epsi library
Epsi_library = '/Users/ncouto/GitHub/EPSILOMETER';
addpath(genpath(Epsi_library));

% Remove archived scripts
rmpath(genpath(fullfile(Epsi_library,'archived')));

% Remove realtime_fctd scripts
rmpath(genpath(fullfile(Epsi_library,'realtime','realtime_epsi')));
