% Find the value closest in one vector (vec) to another given value  (val)

% Bethan Wynne-Cattanach
% 09-01-20

function [closestval,ind] = closest(vec,val)
  ind=NaN(1,length(val));
  for i=1:length(val)
    [~,ind(i)]=min(abs(vec-val(i)),[],'omitnan'); 
  end
  closestval = vec(ind);
end
