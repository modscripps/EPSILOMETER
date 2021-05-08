function [] = check_ctd_time(Meta_Data)

%% Load ctd data
load(fullfile(Meta_Data.CTDpath,['ctd_' Meta_Data.deployment '.mat']));

%% Does ctdtime ever go backwards?
if any(diff(ctd.ctdtime)<0)
 
    dt = mode(diff(ctd.ctdtime));
    n = numel(ctd.ctdtime);
    last_time = n*dt + ctd.ctdtime(1) - dt;
    new_ctdtime = reshape(ctd.ctdtime(1):dt:last_time,n,1);
    

%% Plot figure to check
figure
subplot(2,1,1)
plot(ctd.ctdtime)
title('old ctdtime')

subplot(2,1,2)
plot(new_ctdtime)
title('new ctdtime')

keep = input('Keep new epsitime? (y/n):','s');

%% Update ctdtime
switch keep
    case 'y'
        fprintf('Updating ctd.ctdtime')
        ctd.ctdtime = new_ctdtime;
        save(fullfile(Meta_Data.CTDpath,['ctd_' Meta_Data.deployment '.mat']),'ctd');
    case 'n'
        fprintf('ctd.ctdtime not updated \n')     
end

end