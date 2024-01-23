function [fig] = fullscreenfigure()

fig = figure;

% Set figure size based on screen size
defaultFigWidth = 1680;
defaultFigHeight = 886;
screenSize = get(0,'screensize');
mult = round(min([screenSize(3)/defaultFigWidth,screenSize(4)/defaultFigHeight]),2);
set(fig,'Units','pixels',...
    'Position',[1 1 defaultFigWidth*mult defaultFigHeight*mult],...
    'PaperPosition',[1 1 defaultFigWidth*mult defaultFigHeight*mult]);