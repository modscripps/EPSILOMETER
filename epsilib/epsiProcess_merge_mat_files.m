function MAT = Epsi_MergeMatFiles(MAT1, MAT2)
% MAT = Epsi_MergeMatFiles(MAT1, MAT2)
%
% Nicole Couto edited from FastCTD_MergeFCTD.m
% May 2021
% -------------------------------------------------------
%  MAT = Epsi_MergeMatFiles(MAT1, MAT2)
%   Merges two FCTD structures together keeping all data
%   This also assumes that both FCTD1 and FCTD2 has the same structure in
%   the organization. If not, then an error will be produced and will not
%   be warranted by the code written
%
%   FCTD1 or FCTD2 is empty then FCTD takes the form of the other one
%
%  Written by San Nguyen 2011/07/01
% -------------------------------------------------------

if nargin ~= 2
    error('Must pass in two .mat structures to merge');
end

if isempty(MAT1)
    MAT = MAT2;
    return;
end

if isempty(MAT2)
    MAT = MAT1;
    return;
end

fNames = fieldnames(MAT1);

for i = 1:length(fNames)
    if isempty(MAT1.(fNames{i}))
        if isempty(MAT2.(fNames{i}))
            MAT.(fNames{i}) = [];
        else
            MAT.(fNames{i}) = MAT2.(fNames{i});
        end
        
    % if struct, recurse
    elseif isstruct(MAT1.(fNames{i})) 
        MAT.(fNames{i}) = Epsi_MergeMatFiles(MAT1.(fNames{i}),MAT2.(fNames{i}));
        
    % if data is a cell of strings
    elseif iscellstr(MAT1.(fNames{i})) 
        MAT.(fNames{i}) = MAT1.(fNames{i});
        if ischar(MAT2.(fNames{i}))
            MAT.(fNames{i}){end+1} = MAT2.(fNames{i});
        elseif iscellstr(MAT2.(fNames{i}))
            for j = 1:length(MAT2.(fNames{i}))
                MAT.(fNames{i}){end+1} = MAT2.(fNames{i}){j};
            end
        end
        
    % if this is a string or and array of string (of equal length)
    elseif ischar(MAT1.(fNames{i})) 
        if iscolumn(MAT1.(fNames{i})) || isrow(MAT1.(fNames{i}))
            if ischar(MAT2.(fNames{i}))
                MAT.(fNames{i}) = {MAT1.(fNames{i}); MAT2.(fNames{i})};
            elseif iscellstr(MAT2.(fNames{i}))
                MAT.(fNames{i}) = MAT1.(fNames{i});
                for j = 1:length(MAT2.(fNames{i}))
                    MAT.(fNames{i}){end+1} = MAT2.(fNames{i}){j};
                end
            end
        else
            myDims1 = size(MAT1.(fNames{i}));
            myDims2 = size(MAT2.(fNames{i}));
            if myDims1(1) < myDims1(2)
                if myDims1(1) == myDims2(1)
                    MAT.(fNames{i}) = [MAT1.(fNames{i}) MAT2.(fNames{i})];
                elseif myDims1(2) == myDims2(2)
                    MAT.(fNames{i}) = [MAT1.(fNames{i}); MAT2.(fNames{i})];
                else
                    if ischar(MAT2.(fNames{i}))
                        MAT.(fNames{i}) = {MAT1.(fNames{i}); MAT2.(fNames{i})};
                    elseif iscellstr(MAT2.(fNames{i}))
                        MAT.(fNames{i}) = MAT1.(fNames{i});
                        for j = 1:length(MAT2.(fNames{i}))
                            MAT.(fNames{i}){end+1} = MAT2.(fNames{i}){j};
                        end
                    end
                end
            else
                if myDims1(2) == myDims2(2)
                    MAT.(fNames{i}) = [MAT1.(fNames{i}); MAT2.(fNames{i})];
                elseif myDims1(1) == myDims2(1)
                    MAT.(fNames{i}) = [MAT1.(fNames{i}) MAT2.(fNames{i})];
                else
                    if ischar(MAT2.(fNames{i}))
                        MAT.(fNames{i}) = {MAT1.(fNames{i}); MAT2.(fNames{i})};
                    elseif iscellstr(MAT2.(fNames{i}))
                        MAT.(fNames{i}) = MAT1.(fNames{i});
                        for j = 1:length(MAT2.(fNames{i}))
                            MAT.(fNames{i}){end+1} = MAT2.(fNames{i}){j};
                        end
                    end
                end
            end
        end
    
    % if it is an scalar or an matrix of number
    elseif isnumeric(MAT1.(fNames{i}))
        if length(MAT1.(fNames{i})) > 1
            MAT.(fNames{i}) = [MAT1.(fNames{i}); MAT2.(fNames{i})];
        elseif ~isempty(MAT2.(fNames{i}))
                MAT.(fNames{i}) = [MAT1.(fNames{i}); MAT2.(fNames{i})];
        end
    else
        MAT.(fNames{i}) = [];
    end 
end
return;

end

