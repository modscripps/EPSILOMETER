function vnav=epsiProcess_cpt_shear_heading(vnav)
    % vnav.compass(:,2) is y-axis 
    % vnav.compass(:,x) is x-axis 
    % I am assuming the angle theta to the north would be
    % tan(theta)=y/x   
    vnav.shear_heading=atand(vnav.compass(:,2)./vnav.compass(:,1));
end
