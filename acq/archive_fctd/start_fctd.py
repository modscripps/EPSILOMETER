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



#ser_read = serial.Serial('/dev/tty.usbserial-FTXOY2MY',9600)  # open serial port
#print(ser_read.name)         # check which port was really used
time.sleep(.1)
ser_write = serial.Serial('/dev/tty.usbserial-FT4YHKVU',230400)  # open serial port
print(ser_write.name)         # check which port was really used


utc_time=time.time()

str_cmd  = "time.set %i\r\n" % utc_time


time.sleep(.1)
bin_cmd=str_cmd.encode()
time.sleep(.1)

ser_write.write(bin_cmd)
time.sleep(.1)
print(bin_cmd)


str_cmd  = "settings.mission BLT3 EPSI1 \r\n"
bin_cmd=str_cmd.encode()
ser_write.write(bin_cmd)
time.sleep(.1)
print(bin_cmd)

str_cmd  = "som.fctd\r\n"
bin_cmd=str_cmd.encode()
ser_write.write(bin_cmd)
time.sleep(.1)
print(bin_cmd)

#efe.probe probe_id (1 , 2 ,3 ,4) channel name (t1,t2,s1,s2) probe serial number (3digits) probe calibration (0 for temp, Sv for shear)   
str_cmd  = "efe.probe 1 t1 111 00\r\n" 
bin_cmd=str_cmd.encode()
ser_write.write(bin_cmd)
time.sleep(.1)
print(bin_cmd)


#efe.probe probe_id (1 , 2 ,3 ,4) channel name (t1,t2,s1,s2) probe serial number (3digits) probe calibration (0 for temp, Sv for shear)   
str_cmd  = "efe.probe 2 t2 222 00\r\n" 
bin_cmd=str_cmd.encode()
ser_write.write(bin_cmd)
time.sleep(.1)
print(bin_cmd)


#efe.probe probe_id (1 , 2 ,3 ,4) channel name (t1,t2,s1,s2) probe serial number (3digits) probe calibration (0 for temp, Sv for shear)   
str_cmd  = "efe.probe 3 c1 333 40\r\n"
bin_cmd=str_cmd.encode()
ser_write.write(bin_cmd)
time.sleep(.1)
print(bin_cmd)


#efe.probe probe_id (1 , 2 ,3 ,4) channel name (t1,t2,s1,s2) probe serial number (3digits) probe calibration (0 for temp, Sv for shear)   
str_cmd  = "efe.probe 4 f1 444 40\r\n"
bin_cmd=str_cmd.encode()
ser_write.write(bin_cmd)
time.sleep(.1)
print(bin_cmd)

str_cmd  = "settings.stream\r\n"
bin_cmd=str_cmd.encode()
ser_write.write(bin_cmd)
time.sleep(.1)
print(bin_cmd)



str_cmd  = "som.start\r\n" 
bin_cmd=str_cmd.encode()
ser_write.write(bin_cmd)
time.sleep(.1)
print(bin_cmd)


    
    
    
