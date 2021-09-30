function [Meta_Data] = calibrate_dTdV(Meta_Data);

% A new possible method for finding dTdV for a deployment. Find dTdV for as
% many 50-second overlapping scans as there are in a profile. Pick the best
% one as the median of the top 10% of dTdV values. Plot the 2D PDF to
% visualize the choice.
%
% Nicole Couto | September 2021


fileList = dir(fullfile(Meta_Data.paths.profiles,'Profile*'));
for iFile=1:length(fileList)
        load(fullfile(Meta_Data.paths.profiles,fileList(iFile).name));
        md = mod_epsi_temperature_spectra_v4(Profile.Meta_Data,Profile,0,0);
        dTdV_profile_1{Profile.profNum} = md.AFE.t1.cal_profile(:);
        dTdV_profile_2{Profile.profNum} = md.AFE.t2.cal_profile(:);
        rangeT_profile{Profile.profNum} = md.AFE.t1.ctd_Tmax(:)-md.AFE.t1.ctd_Tmin(:);
        pr_profile{Profile.profNum} = md.AFE.t1.pr(:);
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

[N1,XEDGES,YEDGES] = histcounts2(RANGETall,DTDVall_1,0:0.005:0.15,0.5:1:55,'normalization','pdf');
[N2,XEDGES,YEDGES] = histcounts2(RANGETall,DTDVall_2,0:0.005:0.15,0.5:1:55,'normalization','pdf');
XBINS = nanmean([XEDGES(1:end-1);XEDGES(2:end)]);
YBINS = nanmean([YEDGES(1:end-1);YEDGES(2:end)]);

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

[X,Y] = meshgrid(XBINS,YBINS);
X = X.';
Y = Y.';

dTdV(2) = nanmedian(Y(IDX2_sorted(1:top10_2)));
dTdV(1) = nanmedian(Y(IDX1_sorted(1:top10_1)));

Meta_Data.AFE.t1.cal = dTdV(1);
Meta_Data.AFE.t2.cal = dTdV(2);

save(fullfile(Meta_Data.paths.data,'Meta_Data.mat'),'Meta_Data');

%% Plot t1 dTdV
figure

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
pcolor(XBINS,YBINS,N1.');
shading flat
xlabel('T range')
ylabel('dTdV - t1')
cb(3) = colorbar;
cb(3).Label.String = '% of total observations';

%% Plot t2 dTdV
figure

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
pcolor(XBINS,YBINS,N2.');
shading flat
xlabel('T range')
ylabel('dTdV - t2')
cb(3) = colorbar;
cb(3).Label.String = '% of total observations';