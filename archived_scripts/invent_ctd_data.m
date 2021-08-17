function [] = invent_ctd_data(Meta_Data,pressureOnly)

% [] = invent_ctd_data(Meta_Data,includePressure)
%
% When doing epsi bench tests, we sometimes want to be able to use the full
% EPSILOMETER library to process turbulence data. To do this requires
% pressure, temperature, fall speed data, etc. If we don't have this, we
% make it up here. 
% 
% The number of meters epsi falls in a typical profile is approx 1/35 the
% number of ctd samples and approx. 1/700 the number of epsi samples.

if pressureOnly
    
    % If you're only inventing pressure, that means you have a ctd
    % timeseries from the bench
    % Load ctd data to see how long the timeseries is
    load(fullfile(Meta_Data.CTDpath,['ctd_' Meta_Data.deployment '.mat']),'ctd');
    L = length(ctd.ctdtime);
    pRange = ceil(L/35);
    P = reshape(linspace(0,pRange,L),L,1);
    
    ctd.P = P;
    
    % Also invent dPdt. Turbulence values won't be calculated if w<0.2.
    oneArray = ones(size(P,1),size(P,2));
    ctd.dPdt = 0.5.*oneArray;
    
    save(fullfile(Meta_Data.CTDpath,['ctd_' Meta_Data.deployment '.mat']),'ctd');

elseif ~pressureOnly
    
    % If you're not only inventing pressure, you don't have any ctd data at
    % all
    
    % Load epsi data to see how long the timeseries is
    load(fullfile(Meta_Data.Epsipath,['epsi_' Meta_Data.deployment '.mat']));
    L = length(epsi.epsitime);
    pRange = ceil(L/700);
    P = reshape(linspace(0,pRange,L),L,1);
    
    c3515 = 42.914;
    ctd.P = P;
    oneArray = ones(size(P,1),size(P,2));
    ctd.dzdt = 0.5.*oneArray;
    ctd.T = 12.*oneArray;
    ctd.C = 3.857.*oneArray;
    ctd.S    = sw_salt(ctd.C*10./c3515,ctd.T,ctd.P);
    ctd.sig  = sw_pden(ctd.S,ctd.T,ctd.P,0);
    ctd.ctdtime = reshape(linspace(nanmin(epsi.epsitime),nanmax(epsi.epsitime),L),L,1);
    ctd.ctddnum = reshape(linspace(nanmin(epsi.epsitime),nanmax(epsi.epsitime),L),L,1);
    save(fullfile(Meta_Data.CTDpath,['ctd_' Meta_Data.deployment '.mat']),'ctd') 
    
end
            