function [] = figureStamp(scriptName,varargin)
%
% Usage: figureStamp(matlab.desktop.editor.getActiveFilename)
%        figureStamp(getFilename) <-- getFilename.m is Bethan's easy way to
%        get the current script name
%
% A script that will stamp the current script name and the
% date at the bottom of a figure.
%
% Nicole Couto | March 2019
% -------------------------------------------------------------------------

if numel(varargin)==0
        strPosition = 'bottomleft';
elseif numel(varargin)==1
        strPosition = varargin{1};
end

% Get the full path of the script
slashes = strfind(scriptName,'/');

% The script name will be the name of the script, the directory it's in,
% and the directory that directory is in (so go back three slashes)
if ~isempty(slashes)
    printName = scriptName(slashes(end-2):end);
else
    printName = scriptName;
end

% Get the date
printDate = date;

printString = strrep([printName '.m  ' printDate],'_','\_');

% Find the length of the string and the width of the axes so you can
% position the text box appropriately
char_to_9pt_font_multiplier = 0.004;
str_width_multiplier = 1.4;
w = length(printString)*char_to_9pt_font_multiplier;

axHeightOffset = 0.01;
axWidthOffset = 0.005;

width = w*str_width_multiplier;
height = 0.03;
switch strPosition
    case 'bottomright'
        left = 1 - width;
        bottom = axHeightOffset;
    case 'bottomleft'
        left = axWidthOffset;
        bottom = axHeightOffset;
    case 'topright'
        left = 1 - width;
        bottom = 1 - axHeightOffset - height;
    case 'topleft'
        left = axWidthOffset;
        bottom = 1 - axHeightOffset - height;
end
pos = [left bottom width height];

a = annotation('textbox',pos,'String',printString,'FontSize',9,...
    'EdgeColor','none','BackgroundColor','w','FitBoxToText','on');
