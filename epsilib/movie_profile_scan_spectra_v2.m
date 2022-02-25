function [] = movie_profile_scan_spectra(Profile)
%function mod_epsi_chi_epsi_checkprofile(MS,Meta_Data,l)
% Feb 2019 ALB. Pass in structure MS and Meta_Data.  This will generate a
% standard movie for cast profNum.

%% Get Profile from Profile_or_profNum
% ----------------------------------------------------------
profNum = Profile.profNum;
Meta_Data = Profile.Meta_Data;

%% Initialize movie
% ----------------------------------------------------------

v = VideoWriter(sprintf('%s_cast%i%s',Meta_Data.deployment,profNum,'.avi'));
v.FrameRate=5;
open(v)

for scanNum = 10:10:Profile.nbscan
    
    d = Profile.pr(scanNum);
    
    % Plot the scan data
    [~,~,p5,p6,p7,p8,p9,p10,f] = plot_profile_and_spectra(Profile,d);

    % Grab frame for movie
    pause(.001)
    frame(scanNum)=getframe(gcf);

    
    % Delete scan data, keep things that are constant for profile
    delete(p5(:));
    delete(p6(:));
    delete(p7(:));
    delete(p8(:));
    delete(p9(:));
    delete(p10(:));
    delete(f(:));
   
end %end loop through scans

%% Write and close movie
% ----------------------------------------------------------
writeVideo(v,frame(scanNum))
close(v)

