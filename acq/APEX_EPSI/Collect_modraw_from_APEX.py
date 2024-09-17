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
ser_write = serial.Serial('/dev/tty.usbserial-FT6R60PT',230400,timeout=1)  # open serial port
print(ser_write.name)         # check which port was really used



str_cmd  = "sdio.enable\r\n"
time.sleep(.1)
bin_cmd=str_cmd.encode()
time.sleep(.1)

ser_write.write(bin_cmd)
print(bin_cmd)
readbytes=999
while readbytes>0:
   somline=ser_write.readlines()
   readbytes=ser_write.in_waiting

print(somline)

readbytes=999
while readbytes>0:
   somline=ser_write.readlines()
   readbytes=ser_write.in_waiting

print(somline)

Profile_filename="Profile1_0.modraw"
str_cmd  = "sdio.getdata %s 0\r\n" % Profile_filename
time.sleep(.1)
bin_cmd=str_cmd.encode()
time.sleep(.1)

ser_write.write(bin_cmd)
print(bin_cmd)

# Open a file in write mode
with open(Profile_filename, 'wb') as file:
    try:
        readbytes=999
        while readbytes>0:
            if ser_write.in_waiting > 0:  # Check if there is data available to read
                data = ser_write.readline()  # Read a line of data from the serial port
                file.write(data)  # Decode and write data to file
                file.flush()  # Ensure data is written to file immediately
                print(data)
    except KeyboardInterrupt:
        print("Program interrupted by user.")
    finally:
        ser_write.close()  # Ensure the serial port is properly closed when exiting
        file.close()
    
    
