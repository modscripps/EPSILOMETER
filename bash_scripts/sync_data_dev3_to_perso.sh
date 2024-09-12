#!/bin/bash

while true
do
date +"copy DEV3 to My Drive +%c"
rsync -auq /Volumes/EPSI_PROCESSING/TFO2024/Realtime_RAW/raw/ /Volumes/My\ Drive/DATA/TFO/2024/Realtime_RAW/raw/
rsync -auq /Volumes/EPSI_PROCESSING/TFO2024/Processed/RAW_full_cruise/ /Volumes/My\ Drive/DATA/TFO/2024/Processed/RAW_full_cruise/
sleep 1
done


