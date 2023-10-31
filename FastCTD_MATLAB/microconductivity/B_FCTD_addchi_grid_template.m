%% FCTD Processing (without properly despiking salinity) 
%Bethan Wynne-Cattanach created July 2021
% Updated January 2022 to include inputting timeseries of data rather than
% just MAT files
%% Get the data 
%clear all
%Add FCTD scripts to the path
addpath(genpath('/Users/bethanwynnecattanach/Google Drive/MOD/Processing/FCTD/Microconductivity'))

%Data directory for the FCTD files 
datadir='/Volumes/GoogleDrive-108666672393153374159/Shared drives/MOD data 2021 BLT1/FCTD_EPSI/FCTD_MAT';
    
%What do you want to do? 
%Start from (0) timeseries or (1) MAT files i.e. Merge MAT files?
do_merge=1;
%Calibrate? (0) you can input your own parameters (e.g. if calibration has previously been run)
%  (1) calculates new ones 
do_calibration=1;
%Get chi?
do_chi_calc=1;
%Grid?
do_gridding=1;

%What do you want to save?
%Save merged time series?
save_merged=1;
merged_filename='FCTD_dye_chase_timeseries.mat';
%Save gridded data?
save_gridded=1;
gridded_filename='FCTD_dye_chase_gridded_downcasts_141022.mat';
savedir='/Users/bethanwynnecattanach/Google Drive/MOD/BLT/Data/FCTD/Data/Processed';

%Get chi parameters
chi_param=FCTD_DefaultChiParam;
% Specify the sampling rate of the microconducitivity sensor
chi_param.fs=320; %(320 or 160 are the options accounted for here) 

%% Get files
if do_merge
    %Load the times for each file 
    load(fullfile(datadir,'FastCTD_MATfile_TimeIndex'));

    %%% Specify the time period you're interested in %%%
    times=[datenum(2021,07,01,0,0,0) datenum(2021,07,03,12,0,0)];
    
    [idxFile] = find(FastCTD_MATfile_TimeIndex.timeEnd>times(1) & FastCTD_MATfile_TimeIndex.timeEnd<times(2));
else
    load('FCTD_timeseries.mat')
end

%% Get calibration values for microconductivity 
% This is done for each file individually and then a mean value is
% calculated that can then be used to calculate chi

if do_calibration
    
    if do_merge
    dTdC=nan(1,length(idxFile));
    uC=nan(1,length(idxFile));

    for c=1:length(idxFile)
            idx=idxFile(c);
            fprintf([FastCTD_MATfile_TimeIndex.filenames{idx},'\n'])
            load(fullfile(datadir,FastCTD_MATfile_TimeIndex.filenames{idx}));

            %Reshape the microconductivity to be one long time series
            if chi_param.fs==160
                FCTD.ucon=reshape(FCTD.uConductivity',10*length(FCTD.time),1)/2^16;
            elseif chi_param.fs==320
                FCTD.ucon=reshape(FCTD.uConductivity',20*length(FCTD.time),1)/2^24; 
            else 
                error('Sampling frequency should be 160 or 320Hz')
            end

            %Make the higher resolution time vector for microconductivity
            FCTD.microtime=linspace(FCTD.time(1),FCTD.time(end),chi_param.fs./chi_param.fslow*length(FCTD.time))';

            %Calculate the parameters
            [dTdC(c),uC(c)]=FCTD_ucon_cal(FCTD,chi_param);
            fprintf('dTdC=%.1f, uC=%.1f\n',dTdC(c),uC(c))
    end
    
    %get the average coeffs
    chi_param.gain=nanmean(uC);
    chi_param.dTdC=nanmean(dTdC);
    
    else
        
         %Reshape the microconductivity to be one long time series
            if chi_param.fs==160
                FCTD.ucon=reshape(FCTD.uConductivity',10*length(FCTD.time),1)/2^16;
            elseif chi_param.fs==320
                FCTD.ucon=reshape(FCTD.uConductivity',20*length(FCTD.time),1)/2^24; 
            else 
                error('Sampling frequency should be 160 or 320Hz')
            end

            %Make the higher resolution time vector for microconductivity
            FCTD.microtime=linspace(FCTD.time(1),FCTD.time(end),chi_param.fs./chi_param.fslow*length(FCTD.time))';

            %Calculate the parameters
            [dTdC,uC]=FCTD_ucon_cal(FCTD,chi_param);
            fprintf('dTdC=%.1f, uC=%.1f\n',dTdC(c),uC(c))
    
            chi_param.gain=uC;
            chi_param.dTdC=dTdC;
    end
end

%% Calculate chi 
% This calculates chi for each file and then merges the files together 

if do_chi_calc
    if ~do_calibration
        chi_param.gain=input('gain: ');
        chi_param.dTdC=input('dTdC: ');
    end
    
    if do_merge
    %Load the first file and add chi to initialise the larger merged structure
        load(fullfile(datadir,FastCTD_MATfile_TimeIndex.filenames{idxFile(1)}));
        FCTDall=add_chi_microMHA_v2(FCTD,chi_param);
        for c=2:length(idxFile)
                idx=idxFile(c);
                fprintf([FastCTD_MATfile_TimeIndex.filenames{idx},'\n'])
                load(fullfile(datadir,FastCTD_MATfile_TimeIndex.filenames{idx}));

                 %Reshape the microconductivity to be one long time series
                if chi_param.fs==160
                    FCTD.ucon=reshape(FCTD.uConductivity',10*length(FCTD.time),1)/2^16;
                elseif chi_param.fs==320
                    FCTD.ucon=reshape(FCTD.uConductivity',20*length(FCTD.time),1)/2^24; 
                else 
                    error('Sampling frequency should be 160 or 320Hz')
                end

                %Make the higher resolution time vector for microconductivity
                FCTD.microtime=linspace(FCTD.time(1),FCTD.time(end),chi_param.fs./chi_param.fslow*length(FCTD.time))';

                %Add chi
                FCTD2=add_chi_microMHA_v2(FCTD,chi_param);

                %Merge the files
               FCTDall=FastCTD_MergeFCTD(FCTDall,FCTD2);
        end
    else 
            if chi_param.fs==160
                FCTD.ucon=reshape(FCTD.uConductivity',10*length(FCTD.time),1)/2^16;
            elseif chi_param.fs==320
                FCTD.ucon=reshape(FCTD.uConductivity',20*length(FCTD.time),1)/2^24; 
            else 
                error('Sampling frequency should be 160 or 320Hz')
            end

            %Make the higher resolution time vector for microconductivity
            FCTD.microtime=linspace(FCTD.time(1),FCTD.time(end),chi_param.fs./chi_param.fslow*length(FCTD.time))';

            FCTDall=add_chi_microMHA_v2(FCTD,chi_param);
    end
    
    if save_merged
            save(fullfile(savedir,merged_filename),'FCTDall','-v7.3')
    end
        
end

%% Grid the data

if do_gridding
    %Variables that will be gridded
    % Can be any of the following: {'pressure','temperature','conductivity',
    %'altdist','fluorometer','longitude','latitude','chi','eps','bottomdepth'}
    % Density and Salinity are calculated automatically
    VarsToGrid ={'pressure','temperature','conductivity', 'longitude','latitude', 'chi'};
    
    %Variables that will be concatenated together
    AllVars={[VarsToGrid 'density' 'salinity' 'time']};
    AllVars=AllVars{1,1};

    %Rename/relocate some variables
    FCTDall.longitude=FCTDall.GPS.longitude;
    FCTDall.latitude=FCTDall.GPS.latitude;

    % Get bottom depth if an alitimeter is being used
%     FCTDall.altdist=FCTDall.altDist;
%     % Calculte bottom depth whenever altimeter reads less than 30 m
%      FCTDall.altdist(FCTDall.altdist>30) = nan;
%      depth = sw_dpth(FCTDall.pressure,FCTDall.latitude);
%      FCTDall.bottomdepth = FCTDall.altdist+depth;
%      % Also nan out any value where depth was less than 1000 m
%      FCTDall.bottomdepth(depth<1000) = nan;
%      %clear FCTDall.altDist;
    
    %Take every 10th measurement of the fluorometer
    %FCTDall.fluorometer=((FCTDall.fluorometer(:,1))/2)+1.25;


%     %27/09/22: Calibrating the fluorometer before gridding
%     %Factory calibration:
%     FCTDall.ppb=((FCTDall.fluorometer+2.5)*2-0.0112)*31.1831;
%     FCTDall.ppb(FCTDall.ppb<0)=NaN;
%     FCTDall.ppb=nanmean(FCTDall.ppb,2);


    %Grid the upcasts and the downcasts separately
    FCTDupgrid=B_FastCTD_GridData(FCTDall,'VarsToGrid',VarsToGrid,'upcast');
    FCTDdowngrid=B_FastCTD_GridData(FCTDall,'VarsToGrid',VarsToGrid,'downcast');

    %Concatenate the up and down casts
     numcasts=length(FCTDupgrid.time)+length(FCTDdowngrid.time); %up casts + down casts
     for i=1:length(AllVars)
            varup=FCTDupgrid.(AllVars{i});
            vardown=FCTDdowngrid.(AllVars{i});
            
            % Make all varup bottomdepth = nan. We don't have altimeter
            % data to measure the bottom depth when we're going up
            if strcmp(AllVars{i},'bottomdepth')
               varup =  nan.*FCTDupgrid.(AllVars{i});
            end
            
            
            time_temp=[FCTDdowngrid.time FCTDupgrid.time];
            [~,sort_ind]=sort(time_temp);
            
            var_temp=[vardown varup];
            FCTDgrid.(AllVars{i})=var_temp(:,sort_ind);
     end 
    FCTDgrid.depth=FCTDupgrid.depth;
    FCTDgrid.latitude=nanmean(FCTDgrid.latitude);
    FCTDgrid.longitude=nanmean(FCTDgrid.longitude);
    %FCTDgrid.bottomdepth = nanmean(FCTDgrid.bottomdepth);


    % save the data 
    if save_gridded
        save(fullfile(savedir,gridded_filename),'FCTDgrid','-v7.3')
    end
end
