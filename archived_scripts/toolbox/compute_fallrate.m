function Profile = compute_fallrate(Profile)

%% compute fall rate with pressure
dNH = 2; % for center-difn fall rate, compute using points +/-dNH away
%px = Profile.P(:); tx = Profile.rbrtime(:)*86400;
px = Profile.P(:); tx = Profile.time(:);
fp = dNH+1; lp = length(px)-dNH;
dz = ( px(fp+dNH:lp+dNH) - px(fp-dNH:lp-dNH) ) * 100;
dt = ( tx(fp+dNH:lp+dNH) - tx(fp-dNH:lp-dNH) );
for i = dNH-1:-1:1 % use fewer points near ends
    x(1)=( px(2*i+1)-px(1) )*100; x(2)=( px(end)-px(end-2*i) )*100;
    y(1)=( tx(2*i+1)-tx(1) );  y(2)=( tx(end)-tx(end-2*i) );
    dz = [x(1); dz; x(2)];  dt = [y(1); dt; y(2)];
end
% fwd/bwd difn at ends
dz = [ (px(2)-px(1))*100; dz; (px(end)-px(end-1))*100 ];
dt = [ tx(2)-tx(1); dt; tx(end)-tx(end-1) ];
Profile.w = dz./dt; % FALL RATES (down>0) aligned with pr_scan,time
Profile.dwdt = diff(Profile.w)./diff(tx); 
Profile.dwdt = [Profile.dwdt(1); Profile.dwdt]; % bwd difn accel
