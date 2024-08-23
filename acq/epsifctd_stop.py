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

#tty.usbserial-FTFC08A8

## THIS IS FOR EPSI PC3 McLane

time.sleep(.1)
ser_write = serial.Serial('/dev/tty.usbserial-FTU7BXC3',230400)  # open serial port
print("##connect to: \r\n")
print(ser_write.name)         # check which port was really used


utc_time=time.time()

print("##stop epsi \r\n")

str_cmd  = "som.stop\r\n"
bin_cmd=str_cmd.encode()
ser_write.write(bin_cmd)
time.sleep(.3)
print(bin_cmd)


    
    
    
