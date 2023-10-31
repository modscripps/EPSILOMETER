function X = naninterp(X)
% Interpolate over NaNs
% See INTERP1 for more info
X(isnan(X)) = interp1(find(~isnan(X)), X(~isnan(X)), find(isnan(X)));

if isnan(X(1))==1
    X(1)=X(2);
end
if isnan(X(end))==1
    X(end)=X(end-1);
end

X=interp_sort(X); % make sure values are unique

return