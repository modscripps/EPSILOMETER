#!/usr/bin/env python

"""
DWF Python Epsilometer
Author: Aoife Henry
Revision: 10/17/2013
Inputs: Probe Serial No, Capacitance, Static Head Pressure, Point Skip, Threshold, Dial Indicator Throw
Outputs: Graph of voltage signal output by probe. Calibration Factor,
Requires:
Python 2.7, numpy, matplotlib
python-dateutil, pyparsing

Hard-Coded Parameters: lowPassFrequency, order, widths, # data points to trim from extreme ends of data set, # peak-to-peaks to trim from extreme ends of data set
User-Set Parameters: nSamples, probeSerialNo, capacitanceF capacitanceS, staticHeadPressure, dialIndicatorThrow

"""

# import libraries
from ctypes import *
from dwfconstants import *
import math
import time
import matplotlib.pyplot as plt
import sys
import time
import scipy.signal
import scipy as scp
import numpy as np
import csv
import os

# load dwf libraries
if sys.platform.startswith("win"):
	dwf = cdll.dwf
elif sys.platform.startswith("darwin"):
	dwf = cdll.LoadLibrary("/Library/Frameworks/dwf.framework/dwf")
else:
	dwf = cdll.LoadLibrary("libdwf.so")

# print DWF version
version = create_string_buffer(16)
dwf.FDwfGetVersion(version)
print( "DWF Version: " + str(version.value) )

# raad and check parameters
sv = 27.6473647047 # calculated calibration factor
date = time.strftime("%d/%m/%Y")

inputFlag = False

try:
	
	file = open('config.txt', 'r')
	
	# probe version on line 1 of config.txt
	deviceVersion = int(file.readline().split('= ')[1])

	# number of samples on line 2 of config.txt
	nSamples = int(file.readline().split('= ')[1])
	
	# probe serial number = 5201 on line 3 of config.txt
	probeSerialNo = int(file.readline().split('= ')[1])
	
	# internal capacitance of charge amp feedback path = 986 on line 4 of config.txt
	capacitanceF = float(file.readline().split('= ')[1])
	
	# shear probe float capacitance = 1037 on line 5 of config.txt
	capacitanceS = float(file.readline().split('= ')[1])
	
	# Static Head Pressure = 371 on line 6 of config.txt
	staticHeadPressure = float(file.readline().split('= ')[1])
	
	# Dial Indicator Throw = 0.512 on line 7 of config.txt
	dialIndicatorThrow = float(file.readline().split('= ')[1])

	# Sampling Frequency (samples / sec) on line 8 of config.txt
	hzAcq = c_double( float(file.readline().split('= ')[1]) )

	# Sampling amplitude (V) on line 9 of config.txt
	ampAcq = c_double( float(file.readline().split('= ')[1]) )
	
	# Approximate frequency of shear probe signal (cycles/sec) on line 10 of config.txt
	hzSig = c_double( float(file.readline().split('= ')[1]) )

	# Cutoff frequency for lowpass butterworth filter (cycles/sec) on line 11 of config.txt
	# passes frequencies lower than this value. (Hz). For digital filters, normalized from 0 to 1, where 1 is the Nyquist frequency, pi radians/sample. (thus in half-cycles / sample until converted)
	lowPassFrequency = int(file.readline().split('= ')[1])

	# Order for lowpass butterworth filter (cycles/sec) on line 12 of config.txt. Slightly amplifies signal or produces bad coefficients if too high
	order = int(file.readline().split('= ')[1])
	
	# approximate frequency of noise (cycles/sec) on line 13 of config.txt
	hzNoise = int(file.readline().split('= ')[1])

	file.close()

except:
	print( "Parameters input incorrectly." )
	quit()

## create the folder on mod2 server if needed
data_folder=('/Volumes/MOD_dev/Projects/MMP/SHEAR_PROBE/CALIBRATION/%s' % str(probeSerialNo))
try:
    os.makedirs(data_folder)
except OSError:
    pass


if nSamples < 0:
	
	print( "\nThe value for nSamples is out-of-range;" )
	inputFlag = True

if probeSerialNo < 1:

	print( "\nThe value for Probe serial # is out-of-range;" )
	inputFlag = True

if capacitanceS < 1:
	
	print( "\nThe value for capacitanceS is out-of-range;" )
	inputFlag = True

if staticHeadPressure < 1:
	
	print( "\nThe value for staticHeadPressure is out-of-range;" )
	inputFlag = True

if dialIndicatorThrow < 0.1:
	
	print( "\nThe value for dialIndicatorThrow is out-of-range;" )
	inputFlag = True

if hzAcq.value < 0:
	
	print( "\nThe value for hzAcq is out-of-range;" )
	inputFlag = True

if ampAcq.value < 0:
	
	print( "\nThe value for ampAcq is out-of-range;" )
	inputFlag = True

if hzSig.value < 0.1:
	
	print( "\nThe value for dialIndicatorThrow is out-of-range;" )
	inputFlag = True

if lowPassFrequency < hzSig.value * 10:
	
	print( "\nThe value for lowPassFrequency is out-of-range;" )
	inputFlag = True

if order < 1:
	
	print( "\nThe value for order is out-of-range;" )
	inputFlag = True

if inputFlag == True:

	quit()

# declare ctype variables to facilitate passing of pointers

# Interface handle
hdwf = c_int()

status = c_byte()

# is power supply enabled
isEnabled = c_bool()

# usb voltage reading
usbVoltage = c_double()

# usb current reading
usbCurrent = c_double()

# signal discrete data list
dOut = ( c_double * nSamples )()

# available number of samples to read
cAvailable = c_int()

# number of lost samples after the last check
cLost = c_int()

# number of samples that could be corrupt
cCorrupted = c_int()

# Trigger position to set in seconds
secPosition = c_double( 0 )

# Trigger voltage level to set in volts
#voltsLevel = c_double( 0 )
voltsLevel = c_double( 2 )

# flag indicating if samples were lost
fLost = 0

# flag indicating if samples were corrupt
fCorrupted = 0

# data is input to channel one
channel = c_int(0)

# length of recording time in seconds
recordingLength = nSamples / hzAcq.value

# DATA SOURCE 1: read data from sample_data.csv file.
"""
nSamples = 4096
#recordingLength = 1.0
hzSig = c_double(50)#25)
reader = csv.reader(open('sample_data.csv', 'rb'))
dOut = [ float(x[0]) for x in reader ]
reader.close()
"""

# DATA SOURCE 2: generate signal within python
"""
	
def genSine( A, hzSig, hzAcq, dur ):
	''' Generate a sine wave.
	A = sine wave amplitude
	hzAcq = sample rate (Hz, Samples/sec)
	hzSig = sine wave frequency (Hz)
	duration of signal in seconds = number of samples taken per sec. '''
	t = np.arange( dur )
	sinusoid = A * np.sin( 2 * np.pi * t * ( hzSig / hzAcq ) )
	return sinusoid
	
# sine wave amplitude (V)
A = 1.0
	
hzSig.value = 1
	
# max percentage of sine wave amplitude that noise should have
noiseDepth = 0.1
	
# generate sine wave

sinusoid = genSine( A, hzSig.value, hzAcq.value, nSamples )
	
# generate noise based on normal distrbutions. Mean of the distribution = 0. Scale = Standard deviation. dur = Output shape.
noise = A * noiseDepth * np.random.normal( 0, 0.25, nSamples )
	
dOut = sinusoid + noise
"""

# DATA SOURCE 3: read data from external signal with digilent

# try to open device
print( "Opening first device" )
dwf.FDwfDeviceOpen( c_int( -1 ), byref( hdwf ) )

# if device cannot be opened
if hdwf.value == hdwfNone.value:
	szerr = create_string_buffer(512)
	dwf.FDwfGetLastErrorMsg(szerr)
	print( szerr.value )
	print( "failed to open device" )
	quit()

print( "Preparing to read probe test output signal..." )

# set up analog IO channel nodes
# enable positive supply
dwf.FDwfAnalogIOChannelNodeSet(hdwf, c_int(0), c_int(0), c_double(True))
# set voltage to 2.5 V
#dwf.FDwfAnalogIOChannelNodeSet(hdwf, c_int(0), c_int(1), c_double(2.5))
dwf.FDwfAnalogIOChannelNodeSet(hdwf, c_int(0), c_int(1), c_double(3.3))
# enable negative supply
dwf.FDwfAnalogIOChannelNodeSet(hdwf, c_int(1), c_int(0), c_double(True))
# set voltage to -2.5 V
#dwf.FDwfAnalogIOChannelNodeSet(hdwf, c_int(1), c_int(1), c_double(-2.5))
dwf.FDwfAnalogIOChannelNodeSet(hdwf, c_int(1), c_int(1), c_double(0))
# master enable
dwf.FDwfAnalogIOEnableSet(hdwf, c_int(True))

# set up acquisition
dwf.FDwfAnalogInChannelEnableSet( hdwf, c_int( 0 ), c_bool( True ) )
dwf.FDwfAnalogInChannelRangeSet( hdwf, c_int( 0 ), ampAcq )

# perform acquisition for length of time set by FDwfAnalogInRecordLengthSet.
dwf.FDwfAnalogInAcquisitionModeSet( hdwf, acqmodeRecord )
dwf.FDwfAnalogInFrequencySet( hdwf, hzAcq )
dwf.FDwfAnalogInRecordLengthSet( hdwf, c_double( recordingLength ) )

# set trigger at 0 volts, rising edge
dwf.FDwfAnalogInTriggerSourceSet( hdwf, trigsrcDetectorAnalogIn )
dwf.FDwfAnalogInTriggerPositionSet( hdwf, secPosition )
dwf.FDwfAnalogInTriggerTypeSet( hdwf, trigtypeEdge )
dwf.FDwfAnalogInTriggerConditionSet( hdwf, trigcondRisingPositive )
dwf.FDwfAnalogInTriggerLevelSet( hdwf, voltsLevel )

# wait at least 2 seconds for the offset to stabilize
time.sleep( 2 )

# begin acquisition
dwf.FDwfAnalogInConfigure( hdwf, c_int( 0 ), c_int( 1 ) )
print( "   waiting to finish" )

cSamples = 0

while cSamples < nSamples:
	
	# check supply voltage periodically
	if cSamples % (nSamples/100) == 0:
		#fetch analogIO status from device
		if dwf.FDwfAnalogIOStatus(hdwf) == 0:
			print( "No supply voltage." )
			break

		#supply monitor
		dwf.FDwfAnalogIOChannelNodeStatus(hdwf, c_int(2), c_int(0), byref(usbVoltage))
		dwf.FDwfAnalogIOChannelNodeStatus(hdwf, c_int(2), c_int(1), byref(usbCurrent))
		print "USB: " + str(round(usbVoltage.value,3)) + "V\t" + str(round(usbCurrent.value,3)) + "A"

		# in case of over-current condition the supplies are disabled
		dwf.FDwfAnalogIOEnableStatus(hdwf, byref(isEnabled))
		if not isEnabled:
			#re-enable supplies
			print "Restart"
			dwf.FDwfAnalogIOEnableSet(hdwf, c_int(False))
			dwf.FDwfAnalogIOEnableSet(hdwf, c_int(True))

	dwf.FDwfAnalogInStatus( hdwf, c_int( 1 ), byref( status ) )

	# if Acquisition not yet started.
	if cSamples == 0 and ( status == DwfStateConfig or status == DwfStatePrefill or status == DwfStateArmed ) :

		continue

	dwf.FDwfAnalogInStatusRecord( hdwf, byref( cAvailable ), byref( cLost ), byref( cCorrupted ) )

	cSamples += cLost.value

	if cLost.value :
		fLost = 1

	if cCorrupted.value :
		fCorrupted = 1

	if cAvailable.value == 0 :

		continue

	if cSamples+cAvailable.value > nSamples :

		cAvailable = c_int( nSamples - cSamples )

	# get samples
	dwf.FDwfAnalogInStatusData( hdwf, c_int(0), byref( dOut, 8 * cSamples ), cAvailable )
	cSamples += cAvailable.value

print( "Recording finished" )

# check for errors in sampling process
if fLost:
	print( "Samples were lost! Reduce frequency" )

if cCorrupted:
	print( "Samples could be corrupted! Reduce frequency" )

if nSamples != len( dOut ):
	print( "nSamples not the same length as data recorded...problem?" )

#close the device
dwf.FDwfDeviceClose(hdwf)

# show voltage plot before transformations
fig = plt.figure( figsize=(20,10) )
fig.canvas.set_window_title('Shear Probe Calibrator') 
ax1 = plt.subplot(3,3,1)
plt.plot( np.arange( 0, recordingLength, recordingLength / nSamples), dOut )
plt.title( "Raw Calibration Voltages" )
plt.xlabel( "Time (s)", fontsize=16 )
plt.ylabel( "Voltage (V)", fontsize=16 )

# clone signal
alteredDOut = np.array( dOut )

# approximation of wavelet width based on number of samples between two zero-crossings
waveletWidth = 0

# flag indicating if iterator is still between two zeros points
fWidthMeasurement = 1

# change units to mvs
for i in range ( 0, nSamples ):
	alteredDOut[i] *= 1000

# voltage plot after unit change
ax2 = plt.subplot(3,3,2 )
plt.plot( np.arange(0, recordingLength, recordingLength / nSamples), alteredDOut )
plt.title( "Voltages After Unit Change" )
plt.xlabel( "Time (s)", fontsize=16 )
plt.ylabel( "Voltage (mV)", fontsize=16 )

# calculate average of voltages
avg = sum( alteredDOut ) / nSamples

# calculate standard deviation

# sum of squared elements
ss = 0

for i in range( 0, nSamples ):
	# remove trend
	alteredDOut[i] -= avg
	ss += math.pow( alteredDOut[i], 2 )

# standard deviation
std = ( ss / nSamples )**0.5

# calculate theta TODO???
theta = dialIndicatorThrow / 5.73

# voltage plot after detrending
ax3 = plt.subplot(3,3,3, sharex=ax2, sharey=ax2 )
plt.plot( np.arange(0, recordingLength, recordingLength / nSamples), alteredDOut )
plt.title( "Voltages After Detrending" )
plt.xlabel( "Time (s)", fontsize=16 )
plt.ylabel( "Voltage (mV)", fontsize=16 )

# butterworth lowpass filter

# Nyquist frequency is half the sampling rate.
nyq = 0.5 * hzAcq.value

# normalised cutoff frequency in multiples of nyquiest(Hz)
normalCutoff = lowPassFrequency / nyq

# Numerator (b) and denominator (a) polynomials of the lowpass digital butterworth filter.
b, a = scipy.signal.butter(order, normalCutoff, btype='low', analog=False)

# Compute frequency response of digital filter. w = The normalized frequencies at which h was computed (radians/sample). h = The frequency response.
w, h = scipy.signal.freqz( b, a )

# butterworth filter frequency response plot
fig.add_subplot( 3, 3, 5 )
plt.xscale( "log" )

# convert w from rad/sample to Hz, and h to dBs, to plot
plt.plot( hzAcq.value * w / ( 2.0 * np.pi ), 20 * np.log10( abs( h ) ), 'b')

# label the cutoff frequency and the signal frequency
plt.axvline( lowPassFrequency, color='k' )
plt.plot( lowPassFrequency, 0.5 * np.sqrt(2), 'ko' )
plt.axvline( hzSig.value, color='r' )
plt.axvline( hzNoise, color='g' )
plt.plot( hzSig.value, 0,'ro' )

# set the x limits of the current axes from 0 to nyquist frequency = half of the sampling frequency
plt.xlim( 0, nyq )
plt.title("Butterworth Lowpass Filter Frequency Response")
plt.xlabel("Frequency [Hz]", fontsize=16 )
plt.ylabel("Amplitude [dB]", fontsize=16 )
plt.ylim(-200, 0)

# set x-margin and y-margin
plt.margins( 0, 0.1 )

# display gridlines for both major and minor ticks, and for both x-axis and y-axis
plt.grid( which="both", axis="both" )

# filter alteredDOut with butterworth parameters
alteredDOut = scipy.signal.filtfilt( b, a, alteredDOut )

# voltage plot after filtering
ax4 = plt.subplot(3,3,4, sharex=ax2, sharey=ax2 )
plt.plot( np.arange( 0, recordingLength, recordingLength / nSamples ), alteredDOut )
plt.title( "Voltages After Butterworth Filtering" )
plt.xlabel( "Time (s)", fontsize=16 )
plt.ylabel( "Voltage (mV)", fontsize=16 )

# trim signal. reset nSamples, recordingLength accordingly TODO = WINDOW SIZE
#alteredDOut = np.delete( alteredDOut, range(0, nyq/lowPassFrequency) + range( nSamples + nyq/lowPassFrequency, nSamples) )
nSamples = len(alteredDOut)
recordingLength = nSamples / hzAcq.value

# peak width = peak's full width at half maximum. use small peak width for noisy signal (samples)
# array of width sizes to which the wavelet is stretched to before convolving the wavelet with the data. You should choose a range starting with a value slightly smaller than your expected signal width, up to slightly larger. The more values you supply, the slower the calculation but the higher the resolution. Wavelet = Section of wave from 0 Voltage-Level, to maximum amplitude, back to 0 Voltage-Level

# approximation of wavelet width based on number of samples between two zero-crossings
waveletWidth = 0

# flag indicating if iterator is still between two zeros points
fStartWavelet = 0
fHalfWavelet = 0

# find wavelet width
for i in range( 0 , nSamples - 1 ):
	
	if fStartWavelet == 0:
		if abs(alteredDOut[i+1]) < abs(alteredDOut[i]):
			pass
		elif abs(alteredDOut[i+1]) > abs(alteredDOut[i]):
			fStartWavelet = 1
	
	# if the wavelet is yet to start, and the next zero difference is greater than the old zero difference
	#if fStartWavelet == 0 and abs(newZeroDifference) > abs(oldZeroDifference):
	# if halfway point of wavelet has not yet been reached, and the next zero difference is less than the last zero difference, change flag
	elif fHalfWavelet == 0:
		if abs(alteredDOut[i+1]) > abs(alteredDOut[i]):
			pass
		elif abs(alteredDOut[i+1]) < abs(alteredDOut[i]):
			fHalfWavelet = 1

	# else if the halfway point of the wavelet has been reached, and the next sample slope is greater than the last zero difference, change flag
	# while loop is between two zero-crossings, increment waveletWidth
	elif fStartWavelet == 1 and fHalfWavelet == 1 and abs(alteredDOut[i+1]) > abs(alteredDOut[i]):
		break
	
	waveletWidth += 1

widths = np.arange( round(waveletWidth * 0.95) , round(waveletWidth * 1.05 ))

# find peaks
peakIndices = scipy.signal.find_peaks_cwt( alteredDOut, widths )

# flip all values to opposite sign and find valleys
negativeAlteredDOut = []

for i in range( 0, nSamples ):
	negativeAlteredDOut.append( ( -1 * alteredDOut[ i ] ) )

valleyIndices = scipy.signal.find_peaks_cwt( negativeAlteredDOut, widths )

# initial maximum peak-to-peak
ppMax = 0.0

# peak-to_peak list
numPeakPairs = min( len( peakIndices ), len( valleyIndices ) )
pp = [0.0] * numPeakPairs

for i in range( 0, numPeakPairs ):
	pp[i] = alteredDOut[ peakIndices[i] ] + negativeAlteredDOut[ valleyIndices[i] ]
	if pp[ i ] > ppMax:
		ppMax = pp[ i ]

# remove first and last 2 skewed peak-to-peaks TODO remove outliers instead? Due to incomplete cycles?
pp = pp[ 2:-2 ]
numPeakPairs = len( pp )

# number of valid peak-to-peaks
index = 0

# sum of peak-to-peaks
ppSum = 0.0

for i in range( 0, numPeakPairs ):

	# if the peak-to-peak value is greater then 9.5% of the maximum peak-to-peak voltage found
	if pp[i] > ( 0.095 * ppMax ):

		# increment sum of peak-to-peak voltages
		ppSum += pp[i]
		index += 1

# if there are no peak-to-peaks greater than 9.5% of the maximum peak-to-peak
if index == 0:
	# set average peak-to-peak to 0
	ppAv = 0.0

else:
	# otherwise calculate average peak-to-peak voltage
	ppAv = ppSum / index

# peak-to-peak plot
fig.add_subplot( 3, 3, 6 )
plt.plot( np.arange( 2, numPeakPairs + 2 ), pp )
plt.title( 'Peak to Peak Values' )
plt.xlabel( "Peak-to-Peak Index", fontsize=16 )
plt.ylabel( "Voltage (mV)", fontsize=16 )

# calculate sv
svCalc = ( ( capacitanceF / capacitanceS ) * ppAv ) / ( staticHeadPressure * theta )

# calculate sv difference
svDifference = ( ( svCalc - sv ) / sv ) * 100

# write recorded data to text file

file = open(data_folder + "/record_{}.txt".format( str(probeSerialNo) ), "w")
file.write( "Probe Serial Number = {} \n".format(str(probeSerialNo)) )
file.write( "Probe Version = {} \n".format(str(deviceVersion)) )
file.write( "Internal Capacitance of Charge Amp Feedback Path (uF) = {} \n".format(str(capacitanceF)) )
file.write( "Shear Probe Capacitance (uF) = {} \n".format(str(capacitanceS)) )
file.write( "Static Head Pressure = {} \n".format(str(staticHeadPressure)) )
file.write( "Dial Indicator Throw = {} \n".format(str(dialIndicatorThrow)) )
file.write( "Sampling Frequency (Samples/sec) = {} \n".format(str(hzAcq)) )
file.write( "Sampling Amplitude (V) = {} \n".format(str(ampAcq)) )
file.write( "Approximate Shear Probe Oscillation Frequency (cycles/sec) = {} \n".format(str(hzSig)) )
file.write( "Cutoff Frequency for Lowpass Butterworth Filter (cycles/sec) = {} \n".format(str(lowPassFrequency)) )
file.write( "Order of Lowpass Butterworth Filter = {} \n".format(str(order)) )
file.write( "Approximate Shear Probe Oscillation Noise Frequency (cycles/sec) = {} \n".format(str(hzNoise)) )
file.write( "Number of Samples (Samples) = {} \n".format(str(nSamples)) )
file.write( "Date = {} \n".format(str(date)) )
file.write( "Calculated sv = {} \n".format( str( svCalc ) ) )
file.write( "Original sv = {} \n".format( str( sv ) ) )
file.write( "sv Difference = {}% \n".format( str( svDifference ) ) )
file.write( "Average Peak-to-Peak = {} mV\n".format( str( ppAv ) ) )
file.write( "Standard Deviation of Raw Data = {} mV\n\n".format( str( std ) ) )

for v in alteredDOut:
	file.write( "%s\n" % v )

file.close()

# print calculated values to figure
plt.figtext( 0.5, 0.25, "\n\nAverage Peak-to-Peak is {} mV. \nCalculated sv is {}. Original sv is {}. \nPercentage Difference is {}%. Standard Deviation of Raw Data is {}\n".format( str(ppAv), str( svCalc ), str( sv ), str( svDifference ), str( std ) ), horizontalalignment='center', verticalalignment='baseline' )

# set figure layout and show plot
fig.set_tight_layout(True)

plt.savefig(data_folder +'/Probe' + str(probeSerialNo) +'.png')
plt.show()
