#!/bin/bash

while true
do
date +"copy science party to perso +%c"
rsync -auq /Volumes/science_party_share/RR2410/MOD-SIO/05_processed_data/fctd_epsi/ /Volumes/My\ Drive/DATA/TFO/2024/
sleep 1
done


