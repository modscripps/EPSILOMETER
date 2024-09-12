#!/bin/bash

while true
do
date +"copy RAW to DEV3 +%c"
rsync -auq mod_admin@192.168.1.168:/Users/Shared/FCTD_EPSI_DATA/TFO2024/ /Users/Shared/EPSI_PROCESSING/TFO2024/Realtime_RAW/raw
rsync -auq mod_admin@192.168.1.168:/Users/Shared/FCTD_EPSI_DATA/TFO2024/ /Users/Shared/EPSI_PROCESSING/TFO2024/Processed/RAW_full_cruise
sleep 1 
done


