function [] = movie_profile_scan_spectra(Meta_Data,Profile_or_profNum)
%function mod_epsi_chi_epsi_checkprofile(MS,Meta_Data,l)
% Feb 2019 ALB. Pass in structure MS and Meta_Data.  This will generate a
% standard movie for cast profNum.

%% Get Profile from Profile_or_profNum
% ----------------------------------------------------------

if isnumeric(Profile_or_profNum) && ~isstruct(Profile_or_profNum)
    profNum = Profile_or_profNum;
    load(fullfile(Meta_Data.paths.profiles,sprintf('Profile%03.0f',profNum)));
    eval(['Profile = ' sprintf('Profile%03.0f',profNum) ';']);
elseif isstruct(Profile_or_profNum)
    Profile = Profile_or_profNum;
    profNum = Profile.profNum;
end

%% Initialize movie
% ----------------------------------------------------------

v = VideoWriter(sprintf('%s_cast%i%s',Meta_Data.deployment,profNum,'.avi'));
v.FrameRate=5;
open(v)

for scanNum = 1:Profile.nbscan
    
    % Plot the scan data
    [~,~,p5,p6,p7,p8,p9,p10,f] = plot_profile_scan_spectra(Meta_Data,Profile,scanNum);

    % Grab frame for movie
    pause(.001)
    frame=getframe(gcf);
    writeVideo(v,frame)
    
    % Delete scan data, keep things that are constant for profile
    delete(p5(:));
    delete(p6(:));
    delete(p7(:));
    delete(p8(:));
    delete(p9(:));
    delete(p10(:));
    delete(f(:));
    
end %end loop through scans

%% Close movie
% ----------------------------------------------------------

close(v)
