Meta_Data.mission='PISTON2019';
Meta_Data.vehicle_name='PC2';
Meta_Data.deployment='badinwater';

make_Meta_Data;
Meta_Data.PROCESS.recording_mod='SD';
Meta_Data.PROCESS.nb_channels=8;

Meta_Data.MADRE.rev='MADREB.2';
Meta_Data.MADRE.SN='2';
Meta_Data.MAP.rev='MAPC.0';
Meta_Data.MAP.SN='3';
Meta_Data.CALIpath=fullfile('..','CALIBRATION',[Meta_Data.MADRE.rev '_' ...
    Meta_Data.MADRE.SN '-' Meta_Data.MAP.rev '_' Meta_Data.MAP.SN]);

% TODO: It does not work if the data does not have all the channels
% .     I ll change that if needed

Meta_Data.SBEcal=get_CalSBE('~/ARNAUD/SCRIPPS/EPSILOMETER/SBE49/0131.cal');

Meta_Data.starttime=0;
%% read data
mod_epsi_read_rawfiles(Meta_Data);

%% Calibration
tscan=5;%in seconds
Calibrate_MADREMAP_v2(Meta_Data,tscan);


