matDataDir = '/Users/Shared/EPSI_PROCESSING/Processed/FCTDmat/'; %Directory where FCTD*.mat are stored
rotDataDir = '/Users/Shared/EPSI_PROCESSING/Processed/FCTDrot/'; %Directory where rotation data for each FCTD*.mat file will be stored
%matFiles = dir([matDataDir 'FCTD*.mat']);
matFiles = dir([matDataDir 'EPSI*.mat']);


% if ~exist(rotDataDir,'dir')
%     mkdir(rotDataDir)
% end

file_index_ofset = 1; %index of file in 'matFiles' to start plotting
rot_count_offset = 0; %number of total twists that were on the cable at the end of file_index_offset-1

% -----------------------------------------------------------------------------

% load the index file to indicate what file has been proccessed
if exist([rotDataDir '/matDataFilesLoaded.mat'],'file')
    load([rotDataDir '/matDataFilesLoaded.mat']);
else
    matDataFilesLoaded = repmat(struct('filename',[],'time',[]),0,1);
end

% matDataFilesLoaded{end+1:numel(matFiles)} = cell(1,numel(matFiles)-numel(matDataFilesLoaded));
for i = 1:(numel(matFiles))
    %find out what files are new
    ind = find(strcmpi(matFiles(i).name,{matDataFilesLoaded.filename}));
    if ~isempty(ind)
        %search for time differences
        if(matFiles(i).datenum <= matDataFilesLoaded(ind).time)
            continue;
        end
    end

    %load FCTD data
    %disp(matFiles(i).name);
    load([matDataDir '/' matFiles(i).name]);

    %if the index file is empty, create new index
    if isempty(matDataFilesLoaded)
        matDataFilesLoaded = struct('filename',matFiles(i).name,'time',matFiles(i).datenum);
    else
        if ~isempty(ind)
            matDataFilesLoaded(ind) = struct('filename',matFiles(i).name,'time',matFiles(i).datenum);
        else
            matDataFilesLoaded(end+1) = struct('filename',matFiles(i).name,'time',matFiles(i).datenum);
        end
    end

    % convert mat data file to rotation accumulation file
    if exist('FCTD','var') && isstruct(FCTD) && ~isempty(FCTD.time)
        time = FCTD.time;
        pts = numel(FCTD.time);
%         multiplier = -1;
        multiplier = 1;

        [Phi, Theta, Psi, Rot_mat] = SN_RotateToZAxis([0, 0,1]);
        Rot_Mat = @(p,t,s)[ cos(t)*cos(s), -cos(p)*sin(s) + sin(p)*sin(t)*cos(s),  sin(p)*sin(s) + cos(p)*sin(t)*cos(s);
            cos(t)*sin(s),  cos(p)*cos(s) + sin(p)*sin(t)*sin(s), -sin(p)*cos(s) + cos(p)*sin(t)*sin(s);
            -sin(t),         sin(p)*cos(t),                         cos(p)*cos(t)];
        acce = FCTD.acceleration; %(Rot_mat*(FCTD.acceleration'))';
        acce(:,3) = multiplier*acce(:,3);

        acc_xy_length = sqrt(sum(acce(:,1:2).^2,2));
        acc_length = sqrt(sum(acce.^2,2));

        comp = FCTD.compass;
        gyro = FCTD.gyro;
        acce = FCTD.acceleration;

        comp(:,3) = multiplier*comp(:,3);
        gyro(:,3) = multiplier*gyro(:,3);
        acce(:,3) = multiplier*acce(:,3);
        comp = comp';
        gyro = gyro';
        acce = acce';

        new_comp = NaN(size(comp))';
        new_gyro = NaN(size(comp))';
        new_acce = NaN(size(comp))';
        tic
        for k = 1:pts
            [Phi, Theta, Psi, Rot_mat] = SN_RotateToZAxis(acce(:,k)');
            new_comp(k,:) = Rot_mat*(comp(:,k));
            new_gyro(k,:) = Rot_mat*(gyro(:,k));
            new_acce(k,:) = Rot_mat*(acce(:,k));

        end
        toc
        pressure = FCTD.pressure;
        acceleration = FCTD.acceleration;
        org_gyro = gyro';
        gyro = new_gyro;
        compass = new_comp; %medfilt1(new_comp,5,[],1);
%         compass(end-3:end,:) = new_comp(end-3:end,:);
        COMP_mag = repmat(sqrt(sum(compass(:,1:2).*compass(:,1:2),2)),[1 3]);
        compass = compass./COMP_mag;
        compass = compass(:,1)+1i*compass(:,2);

        dt = diff(time)*24*3600;
        dt = [dt; nanmedian(dt)];

        tot_rot_gyro = cumsum(gyro(:,1).*dt);
        tot_rot_gyro(:,2) = cumsum(gyro(:,2).*dt);
        tot_rot_gyro(:,3) = cumsum(gyro(:,3).*dt);
        tot_rot_acc = phase(compass);
        compass = new_comp;

        save([rotDataDir '/' matFiles(i).name],'compass','tot_rot_gyro','tot_rot_acc','time','org_gyro','gyro','acceleration','pressure');
    end
    clear FCTD;
end
save([rotDataDir 'matDataFilesLoaded.mat'],'matDataFilesLoaded');
%disp('Done');

%%
rotFiles = dir([rotDataDir 'EPSI*.mat']);
% rotFiles = dir([rotDataDir '/deployment_1/FCTD*.mat']);

compass = [];
gyro = [];
org_gyro = [];
tot_rot_gyro = [0, 0, 0];
tot_rot_acc = 0;
time = [];
pressure = [];
acceleration = [];
tic
for i = (1+file_index_ofset):numel(rotFiles)
    %disp(rotFiles(i).name);
    rot = load([rotDataDir '/' rotFiles(i).name]);
    % Make sure there are no nans in the data so you can sum
    rot.tot_rot_gyro(isnan(rot.tot_rot_gyro)) = 0;
    rot.tot_rot_acc(isnan(rot.tot_rot_acc)) = 0;
    % rot = load([rotDataDir '/deployment_1/' rotFiles(i).name]);
    compass = [compass; rot.compass];
    gyro = [gyro; rot.gyro];
    org_gyro = [org_gyro; rot.org_gyro];
    pressure = [pressure; rot.pressure];
    acceleration = [acceleration; rot.acceleration];
    tot_rot_gyro = [tot_rot_gyro; tot_rot_gyro(end,:)+rot.tot_rot_gyro];
    tot_rot_acc = [tot_rot_acc; tot_rot_acc(end,:)+rot.tot_rot_acc];

    time = [time; rot.time];
    clear rot;
end
tot_rot_gyro = tot_rot_gyro(2:end,:);
tot_rot_acc = tot_rot_acc(2:end,:);

toc
disp('Done');

%%
% t_offset = datenum(2019,07,00);
%SN_figure(2,'w',1500,'h',900);
dt = diff(time)*24*3600;
dt = [dt(1); dt];
dp = diff(pressure);
dp = [dp(1); dp];
fall_rate = dp./dt;
% figure(2);
% set(2,'position',[0 0 1500 900])
h = [];
clf;
h(1) = plot(datetime(time,'ConvertFrom','datenum'), -tot_rot_acc/pi/2+rot_count_offset,'linewidth',2,'color','r');
hold on;
h(2) = plot(datetime(time(fall_rate<0.01),'ConvertFrom','datenum'), -tot_rot_acc(fall_rate<0.01)/pi/2+rot_count_offset,'.','linewidth',2,'color','k');

last_value = -tot_rot_acc/pi/2+rot_count_offset;
disp(['Most recent turn count: ' datestr(now),'   ' num2str(round(last_value(end)))])

h(3) = plot(datetime(time,'ConvertFrom','datenum'), tot_rot_gyro(:,3)/pi/2+rot_count_offset,'linewidth',2,'color','b');
h(4) = plot(datetime(time(fall_rate<0.001),'ConvertFrom','datenum'), tot_rot_gyro(fall_rate<0.001,3)/pi/2+rot_count_offset,'.','linewidth',2,'color','c');
set(gca,'XTick',datetime(time(1):12/24:time(end),'ConvertFrom','datenum'))
set(gca,'XTickLabelRotation',45)
hold off;
grid on;
xlabel('time');
ylabel('Number of rotations');
%title('FCTD: Rotation count on current line');
title(['FCTD: Rotation count = ' num2str(round(last_value(end))) '  _ _ _ ']);

hl = legend(h,'Rot by acc [dn]','Rot by acc [up]','Rot by gyro [dn]','Rot by gyro [up]','Location','NorthWest');
set(hl,'Fontsize',20);
set(gca,'Fontsize',20);
ax(1) = gca;
%SN_setTextInterpreter('latex');
% SN_printfig([rotDataDir '/currentRotFig.pdf']);
%
% figure(3);
% set(2,'position',[0 0 1500 900])
% clf;
% plot(datetime(time,'ConvertFrom','datenum'), acceleration,'linewidth',2);
% hold on;
% plot(datetime(time(fall_rate<0.01),'ConvertFrom','datenum'), acceleration(fall_rate<0.01,:),'.','linewidth',2,'color','k');
%
% hold off;
% grid on;
% xlabel('time');
% ylabel('acceleration');
% % title('FCTD: Rotation count on current line');
%
%
% ax(2) = gca;
%
% linkaxes(ax,'x')
