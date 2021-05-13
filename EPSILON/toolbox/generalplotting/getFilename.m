function path=getFilename()
%This is just an easy way to get the filename of the currently open script
%Bethan Wynne-Cattanach March 2021
path = matlab.desktop.editor.getActiveFilename;
end
