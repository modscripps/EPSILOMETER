% I originally created this because I was developing a separate Epsi GUI
% from the FastCTD GUI. The FastCTD GUI uses a bunch of file names that
% call each other and I didn't want to mess with that, but I did want to
% allow the Epsi GUI to plot different variables, so I made two folders
% with the same file names inside. Matlab gets confused and points to the
% wrong files if you have duplicate names like this. In order to switch
% between the two GUI systems to process and plot either Epsi or FastCTD 
% data, I made these 'set_path' scripts to remove one or the other
% directory from the search path.
%
% Long story short, I didn't finish the GUI and never needed this switching
% but I DID like changing the color of my command window for realtime
% visualization vs processing, and for the FastCTD twist counter (which
% does need a different set of paths).
% 
% Nicole Couto | post-BLT3, 2022
% -------------------------------------------------------------------------

% Add Epsi library
startup
Epsi_library = '/Volumes/FCTD Softwares used in BLT 2022/EPSILOMETER_FCTD/';
addpath(genpath(Epsi_library));

% Remove archived scripts
%rmpath(genpath(fullfile(Epsi_library,'archived_scripts')));

% Remove realtime_fctd scripts
%rmpath(genpath(fullfile(Epsi_library,'realtime','realtime_fctd')));

% Set color to yellow
%get properties
cmdWinDoc=com.mathworks.mde.cmdwin.CmdWinDocument.getInstance;
listeners = cmdWinDoc.getDocumentListeners;
%find text area part
jTextArea=listeners(5); %or listeners(3) or listeners (4) depending on matlab
%set colour of command window
jTextArea.setBackground(java.awt.Color.yellow) %for cyan. can also use yellow, pink, etc. and white to turn back 
