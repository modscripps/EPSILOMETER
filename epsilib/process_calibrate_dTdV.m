function [Meta_Data] = process_calibrate_dTdV(Meta_Data);

% A new possible method for finding dTdV for a deployment. Find dTdV for as
% many 50-second overlapping scans as there are in a profile. Pick the best
% one as the median of the top 10% of dTdV values. Plot the 2D PDF to
% visualize the choice.
%
%
% TODO - do this at the .mat file level instead of the Profile level.
%      - the outcome will be almost exactly the same. You might not get as
%      many 50-m segments because some are broken across two files, but
%      you'll be able to calibrate before creating the profiles so you can
%      save a fully-processed profile with turbulence data all in one go.
%      - You might also need to limit to downcasts (or upcasts) to account
%      for a difference due to the relative placement of CTD and t probes
% Nicole Couto | September 2021


fileList = dir(fullfile(Meta_Data.paths.profiles,'Profile*'));
for iFile=1:length(fileList)
    load(fullfile(Meta_Data.paths.profiles,fileList(iFile).name));
    md = mod_epsi_temperature_spectra_v4(Profile.Meta_Data,Profile,0,0);
    if isfield(md,'AFE')
        field_name = 'AFE';
    elseif isfield(md,'epsi')
        field_name = 'epsi';
    end
    dTdV_profile_1{Profile.profNum} = md.(field_name).t1.cal_profile(:);
    dTdV_profile_2{Profile.profNum} = md.(field_name).t2.cal_profile(:);
    rangeT_profile{Profile.profNum} = md.(field_name).t1.ctd_Tmax(:)-md.(field_name).t1.ctd_Tmin(:);
    pr_profile{Profile.profNum} = md.(field_name).t1.pr(:);
end
PROFall= []; PRall=[]; DTDVall_1=[]; DTDVall_2=[]; RANGETall=[];
for ii=1:length(fileList)
    Prof_profile = repmat(ii,length(pr_profile{ii}),1);
    
    DTDVall_1 = [DTDVall_1;dTdV_profile_1{ii}(:)];
    DTDVall_2 = [DTDVall_2;dTdV_profile_2{ii}(:)];
    PROFall = [PROFall;Prof_profile(:)];
    PRall = [PRall;pr_profile{ii}(:)];
    RANGETall = [RANGETall;rangeT_profile{ii}(:)];
end

[~,XEDGES,~] = histcounts(RANGETall(RANGETall<=0.5));
YEDGES = 0.05:1:55;
[N1,XEDGES,YEDGES] = histcounts2(RANGETall(RANGETall<=0.5),DTDVall_1(RANGETall<=0.5),XEDGES,YEDGES,'normalization','pdf');
[N2,XEDGES,YEDGES] = histcounts2(RANGETall(RANGETall<=0.5),DTDVall_2(RANGETall<=0.5),XEDGES,YEDGES,'normalization','pdf');


% ----

lessthan = [0.5:-0.01:0];
for iL=1:length(lessthan)
    med(iL) = nanmedian(DTDVall_1(RANGETall<lessthan(iL)));
    standev(iL) = nanstd(DTDVall_1(RANGETall<lessthan(iL)));
    num(iL) = sum(RANGETall<lessthan(iL));
end
% Sort by where stepping up changes the number of observations the most and
% changes the median dTdV value the least.
[sorted_dnum,rank1] = sort(abs(diff(num)));
[sorted_dmed,rank2] = sort(abs(diff(med)),'descend');
rank1(isnan(sorted_dnum)) = nan;
rank2(isnan(sorted_dnum)) = nan;
rank1 = [nan,rank1];
rank2 = [nan,rank2];
rank = [rank1+rank2];

figure('units','inches','position',[0 0 10 13.1])
subplot(4,1,1)
    plot(lessthan,med,'o-r');
    ylabel('median dTdV')
    title('dTdV - t1')
subplot(4,1,2)
    plot(lessthan,standev,'o-k');
    ylabel('std dTdV')
subplot(4,1,3)
    bar(lessthan,num);
    ylabel('# obs.')
subplot(4,1,4)
    plot(lessthan,[0,abs(diff(num))],'^-k');
    ylabel('change in # obs.')
    xlabel('maximum temperature range')
export_fig figs/dTdV_1_1d -png -r150 -nocrop

% ----

clear med standev num
lessthan = [0.5:-0.01:0];
for iL=1:length(lessthan)
    med(iL) = nanmedian(DTDVall_2(RANGETall<lessthan(iL)));
    standev(iL) = nanstd(DTDVall_2(RANGETall<lessthan(iL)));
    num(iL) = sum(RANGETall<lessthan(iL));
end

figure('units','inches','position',[0 0 10 13.1])
subplot(4,1,1)
plot(lessthan,med,'o-r');
ylabel('median dTdV')
title('dTdV - t2')
subplot(4,1,2)
plot(lessthan,standev,'o-k');
ylabel('std dTdV')
subplot(4,1,3)
bar(lessthan,num);
ylabel('# obs.')
subplot(4,1,4)
plot(lessthan,[0,abs(diff(num))],'^-');
ylabel('change in # obs.')
xlabel('maximum temperature range')
export_fig figs/dTdV_2_1d -png -r150 -nocrop

% ----

maxX   = nanmedian(RANGETall)+2*std(RANGETall);
nX     = round(sqrt(length(PRall)/2));
XEDGES = linspace(0,maxX,nX);
indX   = RANGETall<=maxX;

maxY1 = nanmedian(DTDVall_1)+2*nanstd(DTDVall_1);
minY1 = max([0, nanmedian(DTDVall_1)-2*nanstd(DTDVall_1)]);
maxY2 = nanmedian(DTDVall_2)+2*nanstd(DTDVall_2);
minY2 = max([0, nanmedian(DTDVall_2)-2*nanstd(DTDVall_2)]);
nY   = round(sqrt(length(PRall)*2));
YEDGES1 = linspace(minY1,maxY1,nY);
YEDGES2 = linspace(minY2,maxY2,nY);
indY1 = DTDVall_1>=minY1 & DTDVall_1<=maxY1;
indY2 = DTDVall_2>=minY2 & DTDVall_2<=maxY2;

N1 = histcounts2(RANGETall(indX & indY1),DTDVall_1(indX & indY1),XEDGES,YEDGES1,'normalization','probability');
N2 = histcounts2(RANGETall(indX & indY2),DTDVall_2(indX & indY2),XEDGES,YEDGES2,'normalization','probability');
XBINS = nanmean([XEDGES(1:end-1);XEDGES(2:end)]);
YBINS1 = nanmean([YEDGES1(1:end-1);YEDGES1(2:end)]);
YBINS2 = nanmean([YEDGES2(1:end-1);YEDGES2(2:end)]);

% Find the average values of the top 10% of values
IDX = reshape(1:length(N1(:)),size(N1,1),size(N1,2));
[N1_sorted,i_sort] = sort(N1(:),'descend');
IDX1_sorted = IDX(i_sort);
[N2_sorted,i_sort] = sort(N2(:),'descend');
IDX2_sorted = IDX(i_sort);

N1_cumsum = cumsum(N1_sorted);
N2_cumsum = cumsum(N2_sorted);
[~,top10_1] = min(abs(N1_cumsum-10));
[~,top10_2] = min(abs(N2_cumsum-10));

[X,Y1] = meshgrid(XBINS,YBINS1);
X = X.';
Y1 = Y1.';
[~,Y2] = meshgrid(XBINS,YBINS2);

dTdV(1) = nanmedian(Y1(IDX1_sorted(1:top10_1)));
dTdV(2) = nanmedian(Y2(IDX2_sorted(1:top10_2)));

if isfield(Meta_Data,'AFE')
    Meta_Data.(field_name).t1.cal = dTdV(1);
    Meta_Data.(field_name).t2.cal = dTdV(2);
elseif isfield(Meta_Data,'epsi')
    Meta_Data.epsi.t1.cal = dTdV(1);
    Meta_Data.epsi.t2.cal = dTdV(2);
end

save(fullfile(Meta_Data.paths.data,'Meta_Data.mat'),'Meta_Data');

%% Plot t1 dTdV
figure('units','inches','position',[0 0 10 13.1])

subplot(4,1,1)
scatter(PROFall,PRall,60,DTDVall_1,'filled')
set(gca,'ydir','reverse')
xlabel('Profile #')
ylabel('Pressure')
cb(1) = colorbar;
cb(1).Label.String = 'dTdV - t1';
title(sprintf('t1 - dTdV = %2.1f',dTdV(1)))

subplot(4,1,2)
scatter(PROFall,PRall,60,RANGETall,'filled')
set(gca,'ydir','reverse')
xlabel('Profile #')
ylabel('Pressure')
cb(2) = colorbar;
cb(2).Label.String = 'T range';

subplot(4,1,3:4)
pcolor(XBINS,YBINS1,100*N1.');
shading flat
xlabel('T range')
ylabel('dTdV - t1')
cb(3) = colorbar;
cb(3).Label.String = '% of total observations';

export_fig figs/dTdV_1_2d -png -r150 -nocrop

%% Plot t2 dTdV
figure('units','inches','position',[0 0 10 13.1])

subplot(4,1,1)
scatter(PROFall,PRall,60,DTDVall_2,'filled')
set(gca,'ydir','reverse')
xlabel('Profile #')
ylabel('Pressure')
cb(1) = colorbar;
cb(1).Label.String = 'dTdV - t2';
title(sprintf('t2 - dTdV = %2.1f',dTdV(2)))

subplot(4,1,2)
scatter(PROFall,PRall,60,RANGETall,'filled')
set(gca,'ydir','reverse')
xlabel('Profile #')
ylabel('Pressure')
cb(2) = colorbar;
cb(2).Label.String = 'T range';

subplot(4,1,3:4)
pcolor(XBINS,YBINS2,100*N2.');
shading flat
xlabel('T range')
ylabel('dTdV - t2')
cb(3) = colorbar;
cb(3).Label.String = '% of total observations';

export_fig figs/dTdV_2_2d -png -r150 -nocrop