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

# -----------------------------------------------------------------------------------
port_name="tty.usbserial-FT5953J6"
# -----------------------------------------------------------------------------------

#ser_read = serial.Serial('/dev/tty.usbserial-D308SKZ6',9600)  # open serial port
#print(ser_read.name)         # check which port was really used
time.sleep(.1)
print('/dev/%s' % port_name)         # check which port was really used
ser_write = serial.Serial("/dev/%s"% port_name,230400)   # open serial port
print(ser_write.name)         # check which port was really used

str_cmd  = "som.stop\r\n"

time.sleep(.1)
bin_cmd=str_cmd.encode()
time.sleep(.1)

ser_write.write(bin_cmd)
time.sleep(.1)
print(bin_cmd)
