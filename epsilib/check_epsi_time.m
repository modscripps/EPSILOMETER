function [] = check_epsi_time(Meta_Data)

%% Load epsi data
load(fullfile(Meta_Data.Epsipath,['epsi_' Meta_Data.deployment '.mat']));

%% Does epsitime ever go backwards?
if any(diff(epsi.epsitime)<0)
    

    % Plot figure to check
    ax1 = make_check_time_fig1(epsi);
    ax2 = make_check_time_fig2(epsi);
    
    % Find the spots where it goes backwards.
    dt = mode(diff(epsi.epsitime));
    n = numel(epsi.epsitime);
    last_time = n*dt + epsi.epsitime(1) - dt;
    new_epsitime = reshape(epsi.epsitime(1):dt:last_time,n,1);
    
    keep = input('Keep new epsitime? (y/n):','s');
    
    %% Update epsitime
    switch keep
        case 'y'
            fprintf('Updating epsi.epsitime')
            epsi.epsitime = new_epsitime;
            save(fullfile(Meta_Data.Epsipath,['epsi_' Meta_Data.deployment '.mat']),'epsi');
        case 'n'
            fprintf('epsi.epsitime not updated \n')
    end
    
else
    
     ax1 = make_check_time_fig1(epsi);
     ax2 = make_check_time_fig2(epsi);
    
    
end

end

%% Plots
function [ax] = make_check_time_fig1(epsi)   

figure('units','inches','position',[0 0 10 13])

    ax(1) = subplot(3,1,1);
    plot(epsi.epsitime,'.b')
    shade_different_files(epsi,'samplenum');
    title('epsitime (sec. since power on)')
    
    ax(2) = subplot(3,1,2);
    plot(epsi.efe_time_laptop,'.r')
    shade_different_files(epsi,'samplenum');
    title('laptop time (sec. since ???)')
    
    % Plot epsitime - laptop time (just the first difference of every
    % 80-point block)
    ax(3) = subplot(3,1,3);
    y = epsi.efe_time_laptop(~epsi.efe_bad_blocks).*epsi.efe_block_start(~epsi.efe_bad_blocks);
    plot(epsi.epsitime - y,'.g')
    xlabel('sample number')
    shade_different_files(epsi,'samplenum');
    title('epsitime - laptop time (sec.)')
    
    % Link axes
    linkaxes(ax,'x')
    
end

function [ax] = make_check_time_fig2(epsi)   

    figure('units','inches','position',[0 0 10 13])

    ax(1) = subplot(3,1,1);
    plot(epsi.epsitime - epsi.epsitime(1),'.b')
    shade_different_files(epsi,'samplenum');
    title('relative epsitime (sec.)')
    
    ax(2) = subplot(3,1,2);
    plot(epsi.efe_time_laptop - epsi.efe_time_laptop(1),'.r')
    shade_different_files(epsi,'samplenum');
    title('relative laptop time (sec.)')
    
    % Plot epsitime - laptop time (just the first difference of every
    % 80-point block)
    ax(3) = subplot(3,1,3);
    y = epsi.efe_time_laptop(~epsi.efe_bad_blocks).*epsi.efe_block_start(~epsi.efe_bad_blocks);
    plot((epsi.epsitime - epsi.epsitime(1)) - (epsi.efe_time_laptop(~epsi.efe_bad_blocks).*epsi.efe_block_start(~epsi.efe_bad_blocks) -epsi.efe_time_laptop(1)),'.g')
    xlabel('sample number')
    shade_different_files(epsi,'samplenum');
    title('epsitime - laptop time (sec.)')
    
    % Link axes
    linkaxes(ax,'x')
    
end