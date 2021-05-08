function setpp(dx,dy)
%function setpp(dx,dy)
%Set the x and y widths of the current figure in inches using
%paperposition. 
%MHA 4/11

set(gcf,'paperposition',[1 1 dx dy])
