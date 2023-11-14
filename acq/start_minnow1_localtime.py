#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  9 11:29:33 2020

@author: aleboyer
"""

import numpy as np
from threading import Thread
import serial
import time

#tty.usbserial-FTFC08A8

port_name="tty.usbserial-D308SKZ6"
mission_name="SVALBARD"
vehicle_name="minnow1"
s1_SN="323"
s2_SN="322"
t1_SN="256"
t2_SN="304"

#ser_read = serial.Serial('/dev/tty.usbserial-FTXOY2MY',9600)  # open serial port
#print(ser_read.name)         # check which port was really used
time.sleep(.1)
print('/dev/%s' % port_name)         # check which port was really used

ser_write = serial.Serial("/dev/%s"% port_name,230400)   # open serial port
print(ser_write.name)         # check which port was really used


utc_time=time.time()

str_cmd  = "time.set %i\r\n" % utc_time


time.sleep(.1)
bin_cmd=str_cmd.encode()
time.sleep(.1)

ser_write.write(bin_cmd)
time.sleep(.1)
print(bin_cmd)


str_cmd  = "settings.mission %s %s \r\n" %(mission_name,vehicle_name)
bin_cmd=str_cmd.encode()
ser_write.write(bin_cmd)
time.sleep(.1)
print(bin_cmd)

str_cmd  = "som.epsi\r\n"
bin_cmd=str_cmd.encode()
ser_write.write(bin_cmd)
time.sleep(.1)
print(bin_cmd)

str_cmd  = "som.standalone\r\n"
bin_cmd=str_cmd.encode()
ser_write.write(bin_cmd)
time.sleep(3)
print(bin_cmd)


#efe.probe probe_id (1 , 2 ,3 ,4) channel name (t1,t2,s1,s2) probe serial number (3digits) probe calibration (0 for temp, Sv for shear)   
str_cmd  = "efe.probe 1 t1 %s 00\r\n" % t1_SN
bin_cmd=str_cmd.encode()
ser_write.write(bin_cmd)
time.sleep(.1)
print(bin_cmd)


#efe.probe probe_id (1 , 2 ,3 ,4) channel name (t1,t2,s1,s2) probe serial number (3digits) probe calibration (0 for temp, Sv for shear)   
str_cmd  = "efe.probe 2 t2 %s 00\r\n" % t2_SN
bin_cmd=str_cmd.encode()
ser_write.write(bin_cmd)
time.sleep(.1)
print(bin_cmd)


#efe.probe probe_id (1 , 2 ,3 ,4) channel name (t1,t2,s1,s2) probe serial number (3digits) probe calibration (0 for temp, Sv for shear)   
str_cmd  = "efe.probe 3 s1 %s 40\r\n" % s1_SN
bin_cmd=str_cmd.encode()
ser_write.write(bin_cmd)
time.sleep(.1)
print(bin_cmd)


#efe.probe probe_id (1 , 2 ,3 ,4) channel name (t1,t2,s1,s2) probe serial number (3digits) probe calibration (0 for temp, Sv for shear)   
str_cmd  = "efe.probe 4 s2 %s 40\r\n" %s2_SN
bin_cmd=str_cmd.encode()
ser_write.write(bin_cmd)
time.sleep(.1)
print(bin_cmd)


str_cmd  = "som.start\r\n" 
bin_cmd=str_cmd.encode()
ser_write.write(bin_cmd)
time.sleep(.1)
print(bin_cmd)


    
    
    
