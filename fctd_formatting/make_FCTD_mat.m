function [FCTD] =  make_FCTD_mat(matData,FCTDdir,base,cruise_specifics)

use matData %Empty contents of matData structure

% Get CTD data
FCTD.time=ctd.dnum;
FCTD.pressure=ctd.P;
FCTD.temperature=ctd.T;
FCTD.conductivity=ctd.C;

% Get altimeter data
if ~isempty(alt) && isfield(alt,'time_s')
    FCTD.altDist=interp1(alt.dnum,alt.dst,ctd.dnum);
else
    FCTD.altTime=nan(length(ctd.dnum),1);
end

% Add VectorNav data
if ~isempty(vnav) && isfield(vnav,'time_s')
    diff_not_neg = [0;diff(vnav.dnum)]>0;
    keep = ~isnan(vnav.dnum) & ~isinf(vnav.dnum) & diff_not_neg;
    for ix=1:3
        FCTD.compass(:,ix)=interp1(vnav.dnum(keep),vnav.compass(keep,ix),ctd.dnum);
        FCTD.gyro(:,ix)=interp1(vnav.dnum(keep),vnav.gyro(keep,ix),ctd.dnum);
        FCTD.acceleration(:,ix)=interp1(vnav.dnum(keep),vnav.acceleration(keep,ix),ctd.dnum)./9.81;
    end
else
    FCTD.gyro=nan(length(ctd.dnum),3);
    FCTD.acceleration=nan(length(ctd.dnum),3);
    FCTD.compass=nan(length(ctd.dnum),3);
end

% Add GPS data
if ~isempty(gps) && isfield(gps,'gpstime')
    FCTD.GPS.longitude=interp1(gps.gpstime,gps.longitude,ctd.dnum);
    FCTD.GPS.latitude=interp1(gps.gpstime,gps.latitude,ctd.dnum);
else
    FCTD.GPS.longitude=nan(length(ctd.dnum),1);
    FCTD.GPS.latitude=nan(length(ctd.dnum),1);
end

% Extra outputs for specific cruise setups
if isfield(cruise_specifics,'blt_2021');
    % Microconductivity and Fluorometer
    %
    % On BLT 2021, microconductivity sensor was on shear
    % channel 2 of epsi and fluorometer was on shear channel 1. This step interpolates that data to
    % the same time array as the rest of the data, but since it
    % has a 20x faster sampling rate than the SBE (320 Hz vs 16
    % Hz), it actually becomes and N x 20 array - there are 20
    % uConductivity/fluorometer data points for every 1 SBE data point. We
    % also save time_fast as an N x 20 array.
    time_fast = linspace(ctd.dnum(1),ctd.dnum(end),length(ctd.dnum)*20);
    FCTD.time_fast = time_fast(:);

    % Interpolate data that is not nan, not inf, and where time
    % is increasing
    diff_not_neg = [0;diff(epsi.dnum)]>0;
    keep = ~isnan(epsi.dnum) & ~isinf(epsi.dnum) & diff_not_neg;

    if ~isempty(epsi) && isfield(epsi,'s2_count') && ~isempty(ctd)
        FCTD.uConductivity=reshape(interp1(epsi.dnum(keep),double(epsi.s2_count(keep)),time_fast),20,[])';
    else
        FCTD.uConductivity=nan(length(ctd.dnum),20);
        disp(['No uConductivity data ' myASCIIfiles(i).name]);
    end

    if ~isempty(epsi) && isfield(epsi,'s1_volt')  && ~isempty(ctd)
        FCTD.fluorometer=reshape(interp1(epsi.dnum(keep),epsi.s1_volt(keep),time_fast),20,[])';
    else
        FCTD.fluorometer=nan(length(ctd.dnum),20);
        disp(['No fluorometer data ' myASCIIfiles(i).name]);
    end
end
if isfield(cruise_specifics,'blt_2022')
    
    % On BLT2022, we started adding microstructure values to .mat files
    keep = ~isnan(micro.dnum) & ~isinf(micro.dnum);

    FCTD.epsilon1 = interp1(micro.dnum(keep),micro.epsilon(keep,1),ctd.dnum);
    FCTD.epsilon2 = interp1(micro.dnum(keep),micro.epsilon(keep,2),ctd.dnum);
    FCTD.chi1 = interp1(micro.dnum(keep),micro.chi(keep,1),ctd.dnum);
    FCTD.chi2 = interp1(micro.dnum(keep),micro.chi(keep,2),ctd.dnum);

end %end cruise_specifics

% Save FCTD mat files to the new FCTD mat directory FCTDmat
myFCTDMATfile = fullfile(FCTDdir,base);
save(myFCTDMATfile,'FCTD');
fprintf(1,'%s: Wrote  %s%s\n\n',datestr(now,'YY.mm.dd HH:MM:SS'), FCTDdir,myFCTDMATfile);

% Update FCTD .mat time index
FastCTD_UpdateMATFileTimeIndex(FCTDdir,base,FCTD);

end %end make_FCTD_mat
% ---------------------------------


function FastCTD_UpdateMATFileTimeIndex(dirname,filename,FCTD)
if exist([dirname '/FastCTD_MATfile_TimeIndex.mat'],'file')
    load([dirname '/FastCTD_MATfile_TimeIndex.mat']);
    ind = strncmp(filename,FastCTD_MATfile_TimeIndex.filenames,length(filename));
    if sum(ind) ~= 1
        FastCTD_MATfile_TimeIndex.filenames = [FastCTD_MATfile_TimeIndex.filenames; {filename}];
        FastCTD_MATfile_TimeIndex.timeStart = cat(1,FastCTD_MATfile_TimeIndex.timeStart,FCTD.time(1));
        FastCTD_MATfile_TimeIndex.timeEnd = cat(1,FastCTD_MATfile_TimeIndex.timeEnd,FCTD.time(end));
    else
        FastCTD_MATfile_TimeIndex.timeStart(ind) = FCTD.time(1);
        FastCTD_MATfile_TimeIndex.timeEnd(ind) = FCTD.time(end);
    end
else
    FastCTD_MATfile_TimeIndex.filenames = {filename};
    FastCTD_MATfile_TimeIndex.timeStart = FCTD.time(1);
    FastCTD_MATfile_TimeIndex.timeEnd = FCTD.time(end);
end
save([dirname '/FastCTD_MATfile_TimeIndex.mat'],'FastCTD_MATfile_TimeIndex');
end %end FastCTD_UpdateMATFileTimeIndex
% ------------------------------------