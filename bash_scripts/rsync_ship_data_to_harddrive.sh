#!/bin/sh

rsync -auq --delete /Volumes/data_on_memory/ /Volumes/A-DRIVE/TFO2023/04_ship_data


# Save this script on DEV3
#
# To make the file executable
# chmod a+x copy_raw_to_DEV3
#
# To add to the crontab
# crontab -e
#
# Type the following and save with :w
# 0 */6 * * * /Users/Shared/FCTD_EPSI/rsync_ship_data_to_harddrive
# (The asterisks say when and how often to run the file. This means run every 6 hours at the top of the hour)
# crontab.guru - website to help you figure out how to schedule crontabs 