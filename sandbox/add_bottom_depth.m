function [FCTDgrid] = add_bottom_depth(FCTDgrid)

lat = nanmean(FCTDgrid.latitude(:));

% Add bottom depth from altDist, if you have that
if isfield(FCTDgrid,'altDist')
    for c=219%1:length(FCTDgrid.time)
        if ~isnan(lat)
            ctdZ = sw_dpth(FCTDgrid.pressure(:,c),lat);
        else
            ctdZ = FCTDgrid.pressure(:,c);
        end
        hab = FCTDgrid.altDist(:,c);
        hab(hab>35) = nan;
        bottom_depth = hab+ctdZ;
        % Take the average of the deepest 20 measurements (20 seconds). This step
        % is an attempt to get rid of any spurious readings that might have
        % occurred further up in the profile
        not_nan = ~isnan(bottom_depth);
        bottom_depth = bottom_depth(not_nan);
        ctdZ = ctdZ(not_nan);
        if length(ctdZ)>10
            [~,idxDeep] = sort(ctdZ);
            bottom_depth_median = nanmedian(bottom_depth(idxDeep(end-9:end)));
            FCTDgrid.bottom_depth(c) = bottom_depth_median;
        else
            FCTDgrid.bottom_depth(c) = nan;
        end
    end %End loop through profiles
end