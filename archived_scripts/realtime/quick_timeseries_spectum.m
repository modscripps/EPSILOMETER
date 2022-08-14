%% creates axes
close all
ax(1)=subplot(211);
ax(2)=subplot(212);

%% 20 bit noise
bit20noise=(2.5/2^20)^2/160;

%% 16 bit noise
bit16noise=(2.5/2^16)^2/160;

%% get low epsilon data from BLT
rootpath='/Users/MS/science/APEX/matlab/';
load(fullfile(rootpath,'Low_epsilon_s1_volt_freq_spectrum.mat'));

loglog(ax(2),freq,s1_volt_freq/1000,'b','DisplayName','BLT/1000')
hold(ax(2),'on')
loglog(ax(2),freq,s1_volt_freq,'b--','DisplayName','realBLT:1e-10epsi')
loglog(ax(2),freq,(freq.*0)+bit20noise,'k--','DisplayName','20bit')
loglog(ax(2),freq,(freq.*0)+bit16noise,'k-','DisplayName','16bit')

hold(ax(2),'off')

grid(ax(2),'on')

%% get data from previous tests
rootpath='/Users/MS/science/APEX/matlab/';
listfile=dir(fullfile(rootpath,'*.mat'));
if ~isempty(listfile)
    for f=1:length(listfile)
        if ~strcmp(listfile(f).name,'Low_epsilon_s1_volt_freq_spectrum.mat')
        load(fullfile(rootpath,listfile(f).name));
        hold(ax(2),'on')
        loglog(ax(2),freq,s2_volt_freq,'DisplayName',name)
        hold(ax(2),'off')
        end
    end
end


%% compute current data
ec=epsi_class;
ec.f_readData;
ec.epsi=ec.f_getLastEpsi;

[P1,f1]=pwelch(detrend(ec.epsi.s2_volt(end-21000:end-2500)),[],[],4096,320);

answer0=input("Do you want to plot the new data?(y/n)\n",'s');
if strcmp(answer0,'y')
    plot(ax(1),ec.epsi.time_s(end-11000:end-1000),ec.epsi.s2_volt(end-11000:end-1000))
    hold(ax(2),'on')
    loglog(fftshift(f1),fftshift(P1),'k','linewidth',2,'DisplayName','New test');
    set(ax(2),'Xscale','log','Yscale','log')
    hold(ax(2),'off')
end
    legend(ax(2))
    xlabel(ax(2),'Hz','fontsize',20,'fontname','Times New Roman')
    ylabel(ax(2),'V^2/Hz','fontsize',20,'fontname','Times New Roman')
    xlabel(ax(1),'second','fontsize',20,'fontname','Times New Roman')
    ylabel(ax(1),'Volt','fontsize',20,'fontname','Times New Roman')

    set(ax(1),'fontsize',20,'fontname','Times New Roman');
    set(ax(2),'fontsize',20,'fontname','Times New Roman');

    %% save new data
if strcmp(answer0,'y')
    answer1=input("Do you want to save the new data?(y/n)\n",'s');
    if strcmp(answer1,'y')
        answer2=input("Enter a filename for the new data. default:last_test\n",'s');
        freq=fftshift(f1);
        s2_volt_freq=fftshift(P1);
        name=answer2;
        save(fullfile(rootpath,answer2),'freq','s2_volt_freq','name');
    end
end