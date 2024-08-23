#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  9 11:29:33 2020
modified on Sat 07/30/2022
@author: aleboyer
"""

import numpy as np
from threading import Thread
import serial
import time
import os
import subprocess

###### user input ########

port_name    = 'tty.usbserial-FT4YHELW'
mission_name = "TFO2024"
vehicle_name = "mako1"
s1_SN        = "304"
s2_SN        = "320"
t1_SN        = "326"
t2_SN        = "352"
SBE_SN       = "0537"

##########################



## THIS IS FOR EPSIFISH2 PC2

ser_write = serial.Serial('/dev/%s' % port_name,230400)  # open serial port
time.sleep(.1)
print("##connect to: \r\n")
print(ser_write.name)         # check which port was really used


utc_time=time.time()

print("##start configuring epsi \r\n")

print("##step1 time \r\n")

str_cmd  = "time.set %i\r\n" % utc_time
bin_cmd=str_cmd.encode()
ser_write.write(bin_cmd)
print(bin_cmd)
# read echo from epsi
ser_write.readline()
# read $ from epsi
ser_write.read(1)
answer=ser_write.read(1)
# tell the user the status
if answer!=b'$':
    print("!! step 1 epsi is not responding\r\n")
else:
    print("##step 2 mission\r\n")
time.sleep(.3)
ser_write.read_all()

str_cmd  = "settings.mission %s %s \r\n" %(mission_name,vehicle_name)
bin_cmd=str_cmd.encode()
ser_write.write(bin_cmd)
print(bin_cmd)
# read echo from epsi
ser_write.readline()
# read $ from epsi
ser_write.read(1)
answer=ser_write.read(1)
print(answer)
# tell the user the status
if answer!=b'$':
    print("!! step 2 epsi is NOT responding ")
else:
    print("##step 3 configure ADCs for epsi")
    
ser_write.read_all()
time.sleep(.3)

str_cmd  = "som.epsi\r\n"
bin_cmd=str_cmd.encode()
ser_write.write(bin_cmd)
print(bin_cmd)
# read echo from epsi
ser_write.readline()
# read $ from epsi
ser_write.read(1)
answer=ser_write.read(1)
print(answer)
# tell the user the status
if answer!=b'$':
    print("!! step 3 epsi is NOT responding ")
else:
    print("##step 4 configure probes")
    
ser_write.read_all()
time.sleep(.3)



str_cmd  = "som.telemetry\r\n"
bin_cmd=str_cmd.encode()
ser_write.write(bin_cmd)
print(bin_cmd)
# read echo from epsi
ser_write.readline()
# read $ from epsi
ser_write.read(1)
answer=ser_write.read(1)
print(answer)
# tell the user the status
if answer!=b'$':
    print("!! step 3 epsi is NOT responding ")
else:
    print("##step 4 configure probes")
    
ser_write.read_all()
time.sleep(.3)

#efeieepeee.probe probe_id (1 , 2 ,3 ,4) channel name (t1,t2,s1,s2) probe serial number (3digits) probe calibration (0 for temp, Sv for shear)
str_cmd  = "efe.probe 1 t1 %s 00\r\n" % t1_SN
bin_cmd=str_cmd.encode()
ser_write.write(bin_cmd)
time.sleep(.3)
print(bin_cmd)


#efe.probe probe_id (1 , 2 ,3 ,4) channel name (t1,t2,s1,s2) probe serial number (3digits) probe calibration (0 for temp, Sv for shear)
str_cmd  = "efe.probe 2 t2 %s 00\r\n" %t2_SN
bin_cmd=str_cmd.encode()
ser_write.write(bin_cmd)
time.sleep(.3)
print(bin_cmd)


#efe.probe probe_id (1 , 2 ,3 ,4) channel name (t1,t2,s1,s2) probe serial number (3digits) probe calibration (0 for temp, Sv for shear)
str_cmd  = "efe.probe 3 s1 %s 00\r\n" %s1_SN
bin_cmd=str_cmd.encode()
ser_write.write(bin_cmd)
time.sleep(.3)
print(bin_cmd)



#efe.probe probe_id (1 , 2 ,3 ,4) channel name (t1,t2,s1,s2) probe serial number (3digits) probe calibration (0 for temp, Sv for shear)
str_cmd  = "efe.probe 4 s2 %s 00\r\n" %s2_SN
bin_cmd=str_cmd.encode()
ser_write.write(bin_cmd)
time.sleep(.3)
print(bin_cmd)

str_cmd  = "som.start\r\n"
bin_cmd=str_cmd.encode()
ser_write.write(bin_cmd)
time.sleep(.3)
print(bin_cmd)


ser_write.close()
    
print("starting  fctd_epsi")
time.sleep(2)
subprocess.Popen("/Users/aleboyer/ARNAUD/SCRIPPS/EPSILOMETER/acq/fctd_epsi_acq/build/fctd_epsi/Build/Products/Debug/fctd_epsi")
    
