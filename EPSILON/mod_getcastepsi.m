 function [epsiup,epsidown,dataup,datadown] = mod_getcastepsi(data,CTDtime,up,down)

%  input: data
%   epsi data
%         CTDtime
%  timestamp from the ctd
%         up and down
%  index of the up and down indexes from the casts comouted from the ctd
%
%  output: 
%       epsiup and epsidown
%  epsi indexes fron upcasts and downcasts
%       dataup and datadown
%  epsi data splitted in up/down casts
%
%  Created by Arnaud Le Boyer on 7/28/18.

if isfield(data,'index')
    data=rmfield(data,'index');
    data=rmfield(data,'sderror');
    data=rmfield(data,'chsum1');
    data=rmfield(data,'alti');
    data=rmfield(data,'chsumepsi');
end

 
epsidown=cellfun( @(x) find(data.epsitime>=CTDtime(x(1)) & data.epsitime<=CTDtime(x(end))),down,'un',0);
epsiup=cellfun( @(x) find(data.epsitime>=CTDtime(x(1)) & data.epsitime<=CTDtime(x(end))),up,'un',0);

dataup=cellfun(@(x) structfun(@(y) y(x),data,'un',0),epsiup,'un',0);
datadown=cellfun(@(x) structfun(@(y) y(x),data,'un',0),epsidown,'un',0);

 end
