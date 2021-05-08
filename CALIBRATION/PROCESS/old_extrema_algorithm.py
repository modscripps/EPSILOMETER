isea = 10#nSamples / 100 #int( input( "Enter integer Point Skip: " ) ) # number of steps to wait after a threshold peak-to-peak has been found = 10


if isea < 1:
	
	print( "\nThe value for isea is out-of-range;" )
	inputFlag = True


j = i = 0 # i counter indexes data points. j counter indexes extrema
	
	while i < nSamples:
	
	# maxima search
	extreme.append( threshold )
	
	while i < nSamples:
	
	if alteredDOut[i] < threshold:
	i += 1
	continue
	
	elif alteredDOut[i] > extreme[j]:
	extreme[j] = alteredDOut[i]
	i += 1
	
	elif alteredDOut[i] < extreme[j]:
	break
	
	if i > nSamples - 1:
	break
	
	j += 1
	i += isea
	
	# minima search
	extreme.append( ( -1.0 * threshold ) )
	
	while i < nSamples:
	
	#	if the data point is greater than the threshold, skip it and go to next data point
	if alteredDOut[i] > ( -1.0 * threshold ):
	i += 1
	continue
	
	# else if the data point is less than the existing stored value of the min, reassign min and go to next data point
	elif alteredDOut[i] < extreme[j]:
	extreme[j] = alteredDOut[i]
	i += 1
	
	# else if the data point is greater than the existing stored value of the min, the min has been surpassed. Find next min.
	elif alteredDOut[i] > extreme[j]:
	break
	
	if i > nSamples - 1:
	break
	
	j += 1
	i += isea

# extrema search. End goal is to accurately measure the average peak-to-peak voltage. Finds maxes, find mins. pair consecuitives. find average.
extreme = []

# old extrema list population
extrema = []
m = n = 0
for j in range( len( peakIndices ) ):
	if j % 2 == 0:
		extrema.append( negativeAlteredDOut[ valleyIndices[m] ] )
		m += 1
	else:
		extrema.append( alteredDOut[ peakIndices[n] ] )
		n += 1

pp[i] =extrema[ i * 2 ] + extrema[ ( i * 2 ) + 1 ]

# remove trend alteredDOut = scipy.signal.detrend( alteredDOut )

print(alteredDOut[ peakIndices[0] ])
print(alteredDOut[ peakIndices[numPeakPairs-1] ])
print(negativeAlteredDOut[ valleyIndices[0] ])
print(negativeAlteredDOut[ valleyIndices[numPeakPairs-1] ])


fig.add_subplot( 3, 3, 3 )

def normalise( x, MAX_INT16 ):
	"normalise signal betweem twp amplitude points"
	
	# returns number closest to positive infinity
	maxAmp = max( x )
	
	# MAX_INT16 = maximum signed 16-bit integer value
	amp = math.floor( MAX_INT16 / maxAmp )
	
	# make a list of zeros
	norm = np.zeros( len(x) )
	
	# amplify each data point appropriately
	for i in range( len(x) ):
		norm[i] = amp * x[i]

	return norm

dOut = normalise( dOut, MAX_INT16 )

threshold voltage difference required to qualify as peak-to-peak = 50
threshold = 50 #float( input( "Enter float Threshold: " ) )
if threshold < 1.0:
	
	print( "\nThe value for threshold is out-of-range;" )
	inputFlag = True


# number of samples on line 2 of config.txt
nSamples = int(file.readline().split('= ')[1])#50000 # int( input( "Enter integer Probe Serial Number: " ) ) number of samples
	
	# probe serial number = 5201 on line 3 of config.txt
	probeSerialNo = int(file.readline().split('= ')[1])#5201#int( input( "Enter integer Probe Serial Number: " ) )
	
	# internal capacitance of charge amp feedback path = 986 on line 4 of config.txt
	capacitanceF = float(file.readline().split('= ')[1])#986.0#int( input( "Enter float Capacitance of Charge Amplifier Feedback Path: " ) )
	
	# shear probe float capacitance = 1037 on line 5 of config.txt
	capacitanceS = float(file.readline().split('= ')[1])#1037.0#float( input( "Enter float Capacitance of Shear Probe: " ) )
	
	# Static Head Pressure = 371 on line 6 of config.txt
	staticHeadPressure = float(file.readline().split('= ')[1])#371#float( input( "Enter float Static Head Pressure: " ) )
	
	# Dial Indicator Throw = 0.512 on line 7 of config.txt
	dialIndicatorThrow = float(file.readline().split('= ')[1])#0.512 #float( input( "Enter float Dial Indicator Throw: " ) )

# approximation of wavelet width based on number of samples between two zero-crossings
waveletWidth = 0

# flag indicating if iterator is still between two zeros points
fWidthMeasurement = 1
fHalfWavelet = 0

oldDerivative = (alteredDOut[1] - alteredDOut[0])/(recordingLength*nSamples)
newDerivative = 0

# find wavelet width
for i in range( 1 , nSamples - 1 ):
	newDerivative = (alteredDOut[i + 1] - alteredDOut[i])/(recordingLength*nSamples)
	# if halfway point of wavelet has not yet been reached, and the next sample slope is greater than the last sample slope, change flag
	if fHalfWavelet == 0 and abs(newDerivative) > abs(oldDerivative):
		fHalfWavelet = 1
	# else if the halfway point of the wavelet has been reached, and the next sample slope is less than the last sample slope, change flag
	elif fHalfWavelet == 1 and abs(newDerivative) < abs(oldDerivative):
		fWidthMeasurement = 0
	# while loop is between two zero-crossings, increment waveletWidth
	if fWidthMeasurement == 1:
		waveletWidth += 1
