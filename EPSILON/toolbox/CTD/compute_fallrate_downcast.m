function w = compute_fallrate_downcast(Profile)
%% compute fall rate with pressure
px = Profile.P(:); tx = Profile.ctdtime(:)*86400;
w=diff(px)./nanmean(diff(tx));
w=[w;w(end)];
% 
