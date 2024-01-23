function [] = set_window_color(color)
% set_window_color(color)
%
% color = 'yellow','cyan','pink','white','gray','green'

cmdWinDoc=com.mathworks.mde.cmdwin.CmdWinDocument.getInstance;
listeners = cmdWinDoc.getDocumentListeners;
%find text area part
jTextArea=listeners(5); %or listeners(3) or listeners (4) depending on matlab
%set colour of command window
jTextArea.setBackground(java.awt.Color.(color)) %for cyan. can also use yellow, pink, etc. and white to turn back 
