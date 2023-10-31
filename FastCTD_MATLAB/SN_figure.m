function SN_figure(varargin)
%SN_FIGURE is an extension of FIGURE where it also adds the source of
% creation of the figure
%
%   SN_FIGURE, by itself, uses FIGURE to create a new figure window, 
%   and returns its handle.
% 
%   SN_FIGURE(H) makes H the current figure, forces it to become visible,
%   and raises it above all other figures on the screen.  If Figure H
%   does not exist, and H is an integer, a new figure is created with
%   handle H.
%   
%   SN_FIGURE(...,'OPTIONS',OPT_VAL) allows users to add more advance
%   options
%   The OPTIONS are:
%       'MFILENAME' or 'FNAME':  the name of the m-file that execute the
%           figure
%       'MFILELINE' or 'LINE':  the line number from the the m-file that
%           execute the figure
%       'WIDTH' or 'W': width of figure
%       'HEIGHT' or 'H': height of figure
%       'X', 'Y': x or y position of figure, respectively
%
% See also FIGURE, SN_PRINTFIG
%
% Created by San Nguyen 2014 10 26
% Updated by San Nguyen 2014 11 05
% Updated by San Nguyen 2015 09 13 % added more options
%

mFileName = '';
mFileLine = NaN;
width = NaN;
height = NaN;
x = NaN;
y = NaN;
fname = NaN;
persistent argsNameToCheck;
if isempty(argsNameToCheck);
    argsNameToCheck = {...
        'mFileName',... %1
        'mFileLine',... %2
        'line',... %3
        'fName',... %4
        'width',... %5
        'height',... %6
        'w',... %7
        'h',... %8
        'x',... %9
        'y',... %10
        'position',... %11
        'pos',...%12
        'dim',...%13
        'size',... %14
        'name',... %15
        'figure name'... %16
        };
                        
end

index = 1;
n_items = nargin;
otherVarargin = {};
j = 0;
while (n_items > 0)
    argsMatch = strcmpi(varargin{index},argsNameToCheck);
    i = find(argsMatch,1);
    if isempty(i)
        j = j+1;
        otherVarargin{j} = varargin{index};
        index = index +1;
        n_items = n_items-1;
        continue;
    end
    
    switch i
        case {1,4} % mFileName
            if n_items == 1
                error('MATLAB:SN_figure:missingArgs',...
                    'Missing input arguments');
            end
            mFileName = varargin{index+1};
            if isempty(mFileName) || ~ischar(mFileName)
                error('MATLAB:SN_figure:mFileName',...
                    'Please check your mFileName');
            end
            index = index +2;
            n_items = n_items-2;
        case {2,3} % mFileLine
            if n_items == 1
                error('MATLAB:SN_figure:missingArgs',...
                    'Missing input arguments');
            end
            mFileLine = varargin{index+1};
            if ~(isscalar(mFileLine) && isnumeric(mFileLine) && mod(mFileLine,1)==0)
                error('MATLAB:SN_figure:mFileLine',...
                    'Please check your mFileLine');
            end
            index = index +2;
            n_items = n_items-2;
        case {5,7} % width
            if n_items == 1
                error('MATLAB:SN_figure:missingArgs',...
                    'Missing input arguments');
            end
            width = varargin{index+1};
            if ~(isscalar(width) && isnumeric(width) && mod(width,1)==0 && width > 0)
                error('MATLAB:SN_figure:width',...
                    'Please check your width. Must be an interger');
            end
            index = index +2;
            n_items = n_items-2;
        case {6,8} % height
            if n_items == 1
                error('MATLAB:SN_figure:missingArgs',...
                    'Missing input arguments');
            end
            height = varargin{index+1};
            if ~(isscalar(height) && isnumeric(height) && mod(height,1)==0 && height > 0)
                error('MATLAB:SN_figure:height',...
                    'Please check your height. Must be an interger');
            end
            index = index +2;
            n_items = n_items-2;
        case 9 % x
            if n_items == 1
                error('MATLAB:SN_figure:missingArgs',...
                    'Missing input arguments');
            end
            x = varargin{index+1};
            if ~(isscalar(x) && isnumeric(x)  && mod(x,1)==0 && x > 0)
                error('MATLAB:SN_figure:x',...
                    'Please check your x position. Must be an interger');
            end
            index = index +2;
            n_items = n_items-2;
        case 10 % y
            if n_items == 1
                error('MATLAB:SN_figure:missingArgs',...
                    'Missing input arguments');
            end
            y = varargin{index+1};
            if ~(isscalar(y) && isnumeric(y)  && mod(y,1)==0 && y > 0)
                error('MATLAB:SN_figure:y',...
                    'Please check your y position. Must be an interger');
            end
            index = index +2;
            n_items = n_items-2;
        case {11,12} % position
            if n_items == 1
                error('MATLAB:SN_figure:missingArgs',...
                    'Missing input arguments');
            end
            position = varargin{index+1};
            if isvector(position) || numel(position)~=4 || isnumeric(position) || sum(position<=0)
                error('MATLAB:SN_figure:position',...
                    'Please check your position. It must be a vector of 4 elements.');
            end
            x = position(1);
            y = position(2);
            width = position(3);
            height = position(4);
            if ~(isscalar(height) && isnumeric(height) && mod(height,1)==0 && height > 0) ||...
                    ~(isscalar(width) && isnumeric(width) && mod(width,1)==0 && width > 0) ||...
                    ~(isscalar(x) && isnumeric(x)  && mod(x,1)==0 && x > 0) || ...
                    ~(isscalar(y) && isnumeric(y)  && mod(y,1)==0 && y > 0)
                error('MATLAB:SN_figure:position',...
                    'Please check your position input.');
            end
            index = index +2;
            n_items = n_items-2;  
        case {13,14} % size, dimension
            if n_items == 1
                error('MATLAB:SN_figure:missingArgs',...
                    'Missing input arguments');
            end
            dim = varargin{index+1};
            if isvector(dim) || numel(dim)~=2 || isnumeric(dim) || sum(dim<=0)
                error('MATLAB:SN_figure:size',...
                    'Please check your size values. It must be a vector of 2 elements.');
            end
            width = dim(1);
            height = dim(2);
            if ~(isscalar(height) && isnumeric(height) && mod(height,1)==0 && height > 0) ||...
                    ~(isscalar(width) && isnumeric(width) && mod(width,1)==0 && width > 0) ||...
                error('MATLAB:SN_figure:size',...
                    'Please check your size input.');
            end
            index = index +2;
            n_items = n_items-2; 
        case {15,16} % figure name
            if n_items == 1
                error('MATLAB:SN_figure:missingArgs',...
                    'Missing input arguments');
            end
            fname = varargin{index+1};
            if isvector(fname) || ischar(fname)
                error('MATLAB:SN_figure:figureName',...
                    'Figure name must be a string');
            end
            index = index +2;
            n_items = n_items-2; 
    end
end


fig = figure(otherVarargin{:});
[ST,~] = dbstack('-completenames');

for i=1:numel(ST)
    if strncmp(ST(i).name,'SN_figure',9)
        continue;
    else
        tmp = ST(i).name;
        ind = strfind(tmp,'/');
        mFileLine = ST(i).line;
        if isempty(ind)
            mFileName = tmp;
        else
            mFileName = tmp(1:ind(1)-1);
        end
        break;
    end
end

figDate = now;
if ~isempty(mFileName)
    set(fig,...
        'UserData',struct('mfilename',mFileName,'line',mFileLine));
end
pos = get(fig,'position');

if ~isnan(x)
    pos(1) = x;
end

if ~isnan(y)
    pos(2) = y;
end

if ~isnan(width)
    pos(3) = width;
end

if ~isnan(height)
    pos(4) = height;
end

set(fig,'position',pos);

if ~isnan(fname)
    set(fig,'name',fname);
end
end