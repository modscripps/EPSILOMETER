% Add Epsi library
startup
Epsi_library = '/Volumes/FCTD Softwares used in BLT 2022/EPSILOMETER_FCTD/';
addpath(genpath(Epsi_library));

% Remove archived scripts
rmpath(genpath(fullfile(Epsi_library,'archived_scripts')));

% Remove realtime_fctd scripts
rmpath(genpath(fullfile(Epsi_library,'realtime','realtime_fctd')));

% Set color to yellow
%get properties
cmdWinDoc=com.mathworks.mde.cmdwin.CmdWinDocument.getInstance;
listeners = cmdWinDoc.getDocumentListeners;
%find text area part
jTextArea=listeners(5); %or listeners(3) or listeners (4) depending on matlab
%set colour of command window
jTextArea.setBackground(java.awt.Color.yellow) %for cyan. can also use yellow, pink, etc. and white to turn back 
