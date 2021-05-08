def findpeaksb(x,y,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,smoothtype,window,PeakShape,extra,NumTrials,AUTOZERO):

		sizex=size(x)
		if sizex(1)>sizex(2),x=x'
		sizey=size(y)
		if sizey(1)>sizey(2),y=y'
		NumPeaks=1 # Number of models shapes to you to fit each peak (usually 1)
		start=0
		fixedparameters=0 # Used for fixed parameters shapes only (11,16,12,and 17)
		plots=0 # Supress plotting
		if smoothtype>3smoothtype=3
		if smoothtype<1smoothtype=1 
		smoothwidth=round(smoothwidth)
		peakgroup=round(peakgroup)
		if smoothwidth>1,
			d=fastsmooth(deriv(y),smoothwidth,smoothtype)
		else:
			d=y

		n=round(peakgroup/2+1)
		P=[0 0 0 0 0 0 0]
		vectorlength=length(y)
		peak=1
		AmpTest=AmpThreshold
		for j=2*round(smoothwidth/2)-1:length(y)-smoothwidth,
			if sign(d(j)) > sign (d(j+1)), # Detects zero-crossing
					if d(j)-d(j+1) > SlopeThreshold, # if slope of derivative is larger than SlopeThreshold
							if y(j) > AmpTest,  # if height of peak is larger than AmpThreshold
									xx=zeros(size(peakgroup))yy=zeros(size(peakgroup))
									for k=1:peakgroup, # Create sub-group of points near peak
											groupindex=j+k-n+1
											if groupindex<1, groupindex=1
											if groupindex>vectorlength, groupindex=vectorlength
											xx(k)=x(groupindex)yy(k)=y(groupindex)
									
									if peakgroup>3,
											startpoint=round(j-window/2)
											if startpoint<1startpoint=1
											point=round(j+window/2)
											if point>length(x)point=length(x) 
											XXX=x((startpoint):(point))
											YYY=y((startpoint):(point))  
											signal=[XXX',YYY']
											[FitResults,Error]=peakfit(signal,x(j),window,NumPeaks,PeakShape,extra,NumTrials,start,AUTOZERO,fixedparameters,plots)                
											PeakX=FitResults(2)
											PeakY=FitResults(3)
											MeasuredWidth=FitResults(4)
									else:
											# if the peak is too narrow for least-squares technique
											# to work well, just use the max value of y in the
											# sub-group of points near peak.
											PeakY=max(yy)
											pindex=val2ind(yy,PeakY)
											PeakX=xx(pindex(1))
											MeasuredWidth=0
									
									# Construct matrix P. One row for each peak detected,
									# containing the peak number, peak position (x-value) and
									# peak height (y-value). If peak measurements fails and
									# results in NaN, skip this peak
									if isnan(PeakX) || isnan(PeakY) || PeakY<AmpThreshold,
											# Skip any peak that gives a NaN or is less that the
											# AmpThreshold
									else: # Otherwiase count this as a valid peak
											P(peak,:) = [round(peak) PeakX PeakY MeasuredWidth  1.0646.*PeakY*MeasuredWidth Error]
											peak=peak+1 # Move on to next peak
									 # if isnan(PeakX)...
							 # if y(j) > AmpTest...
						# if d(j)-d(j+1) > SlopeThreshold...
				# if sign(d(j)) > sign (d(j+1))...
		# for j=....
		# ----------------------------------------------------------------------
		def [index,closestval]=val2ind(x,val)
		# Returns the index and the value of the element of vector x that is closest to val
		# If more than one element is equally close, returns vectors of indicies and values
		# Tom O'Haver (toh@umd.edu) October 2006
		# Examples: If x=[1 2 4 3 5 9 6 4 5 3 1], then val2ind(x,6)=7 and val2ind(x,5.1)=[5 9]
		# [indices values]=val2ind(x,3.3) returns indices = [4 10] and values = [3 3]
		dif=abs(x-val)
		index=find((dif-min(dif))==0)
		closestval=x(index)

		def d=deriv(a)
		# First derivative of vector using 2-point central difference.
		#  T. C. O'Haver, 1988.
		n=length(a)
		d(1)=a(2)-a(1)
		d(n)=a(n)-a(n-1)
		for j = 2:n-1
		d(j)=(a(j+1)-a(j-1)) ./ 2


		def SmoothY=fastsmooth(Y,w,type,s)
		# fastbsmooth(Y,w,type,s) smooths vector Y with smooth 
		#  of width w. Version 2.0, May 2008.
		# The argument "type" determines the smooth type:
		#   If type=1, rectangular (sliding-average or boxcar) 
		#   If type=2, triangular (2 passes of sliding-average)
		#   If type=3, pseudo-Gaussian (3 passes of sliding-average)
		# The argument "s" controls how the "s" of the signal 
		# (the first w/2 points and the last w/2 points) are handled.
		#   If s=0, the s are zero.  (In this mode the elapsed 
		#     time is indepent of the smooth width). The fastest.
		#   If s=1, the s are smoothed with progressively 
		#     smaller smooths the closer to the . (In this mode the  
		#     elapsed time increases with increasing smooth widths).
		# fastsmooth(Y,w,type) smooths with s=0.
		# fastsmooth(Y,w) smooths with type=1 and s=0.
		# Example:
		# fastsmooth([1 1 1 10 10 10 1 1 1 1],3)= [0 1 4 7 10 7 4 1 1 0]
		# fastsmooth([1 1 1 10 10 10 1 1 1 1],3,1,1)= [1 1 4 7 10 7 4 1 1 1]
		#  T. C. O'Haver, May, 2008.
		if nargin==2, s=0 type=1 
		if nargin==3, s=0 
		switch type
			case 1
				 SmoothY=sa(Y,w,s)
			case 2   
				 SmoothY=sa(sa(Y,w,s),w,s)
			case 3
				 SmoothY=sa(sa(sa(Y,w,s),w,s),w,s)


		def SmoothY=sa(Y,smoothwidth,s)
		w=round(smoothwidth)
		SumPoints=sum(Y(1:w))
		s=zeros(size(Y))
		halfw=round(w/2)
		L=length(Y)
		for k=1:L-w,
		 s(k+halfw-1)=SumPoints
		 SumPoints=SumPoints-Y(k)
		 SumPoints=SumPoints+Y(k+w)

		s(k+halfw)=sum(Y(L-w+1:L))
		SmoothY=s./w
		# Taper the s of the signal if s=1.
		if s==1,
			startpoint=(smoothwidth + 1)/2
			SmoothY(1)=(Y(1)+Y(2))./2
			for k=2:startpoint,
				 SmoothY(k)=mean(Y(1:(2*k-1)))
				 SmoothY(L-k+1)=mean(Y(L-2*k+2:L))
			
			SmoothY(L)=(Y(L)+Y(L-1))./2

		# ----------------------------------------------------------------------
		def [FitResults,GOF,baseline,coeff,BestStart,xi,yi,BootResults]=peakfit(signal,center,window,NumPeaks,peakshape,extra,NumTrials,start,autozero,fixedparameters,plots,bipolar,minwidth,DELTA,clipheight)
		# A command-line peak fitting program for time-series signals.
		# Version 7: March, 2015, adds peak shapes with three unconstrained
		# iterated variables: 30=voigt (variable alpha), 31=ExpGaussian (variable
		# time constant), 32=Pearson (variable shape factor), 34=Gaussian/
		# Lorentzian bl (variable percent).
		# 
		global AA xxx PEAKHEIGHTS FIXEDPARAMETERS AUTOZERO delta BIPOLAR CLIPHEIGHT
		format short g
		format compact
		warning off all
		NumArgOut=nargout
		datasize=size(signal)
		if datasize(1)<datasize(2),signal=signal'
		datasize=size(signal)
		if datasize(2)==1, #  Must be isignal(Y-vector)
			X=1:length(signal) # Create an indepent variable vector
			Y=signal
		else:
			# Must be isignal(DataMatrix)
			X=signal(:,1) # Split matrix argument 
			Y=signal(:,2)

		X=reshape(X,1,length(X)) # Adjust X and Y vector shape to 1 x n (rather than n x 1)
		Y=reshape(Y,1,length(Y))
		# If necessary, flip the data vectors so that X increases
		if X(1)>X(length(X)),
			disp('X-axis flipped.')
			X=fliplr(X)
			Y=fliplr(Y)

		# Isolate desired segment from data set for curve fitting
		if nargin==1 || nargin==2,center=(max(X)-min(X))/2window=max(X)-min(X)
		# Y=Y-min(Y)
		xoffset=0
		n1=val2ind(X,center-window/2)
		n2=val2ind(X,center+window/2)
		if window==0,n1=1n2=length(X)
		xx=X(n1:n2)-xoffset
		yy=Y(n1:n2)
		ShapeString='Gaussian'
		coeff=0
		CLIPHEIGHT=max(Y)
		LOGPLOT=0
		# Define values of any missing arguments
		# (signal,center,window,NumPeaks,peakshape,extra,NumTrials,start,autozero,fixedparameters,plots,bipolar,minwidth,DELTA)
		switch nargin
			case 1
					NumPeaks=1
					peakshape=1
					extra=0
					NumTrials=1
					xx=Xyy=Y
					start=calcstart(xx,NumPeaks,xoffset)
					AUTOZERO=0
					plots=1
					BIPOLAR=0
					MINWIDTH=xx(2)-xx(1)
					delta=1
					CLIPHEIGHT=max(Y)
			case 2
					NumPeaks=1
					peakshape=1
					extra=0
					NumTrials=1
					xx=signalyy=center
					start=calcstart(xx,NumPeaks,xoffset)
					AUTOZERO=0
					plots=1
					BIPOLAR=0
					MINWIDTH=xx(2)-xx(1)
					delta=1
					CLIPHEIGHT=max(Y)
			case 3
					NumPeaks=1
					peakshape=1
					extra=0
					NumTrials=1
					start=calcstart(xx,NumPeaks,xoffset)
					AUTOZERO=0
					FIXEDPARAMETERS=0
					plots=1
					BIPOLAR=0
					MINWIDTH=xx(2)-xx(1)
					delta=1
					CLIPHEIGHT=max(Y)
			case 4 # Numpeaks specified in arguments
					peakshape=1
					extra=0
					NumTrials=1
					start=calcstart(xx,NumPeaks,xoffset)
					AUTOZERO=0
					FIXEDPARAMETERS=0
					plots=1
					BIPOLAR=0
					MINWIDTH=xx(2)-xx(1)
					delta=1
					CLIPHEIGHT=max(Y)
			case 5 # Numpeaks, peakshape specified in arguments
					extra=zeros(1,NumPeaks)
					NumTrials=1
					start=calcstart(xx,NumPeaks,xoffset)
					AUTOZERO=0
					FIXEDPARAMETERS=0
					plots=1
					BIPOLAR=0
					MINWIDTH=zeros(size(peakshape))+(xx(2)-xx(1))
					delta=1
					CLIPHEIGHT=max(Y)
			case 6
					NumTrials=1
					start=calcstart(xx,NumPeaks,xoffset)
					AUTOZERO=0
					FIXEDPARAMETERS=0
					plots=1
					BIPOLAR=0
					MINWIDTH=zeros(size(peakshape))+(xx(2)-xx(1))
					delta=1
			case 7
					start=calcstart(xx,NumPeaks,xoffset)
					AUTOZERO=0
					FIXEDPARAMETERS=0
					plots=1
					BIPOLAR=0
					MINWIDTH=zeros(size(peakshape))+(xx(2)-xx(1))
					delta=1
					CLIPHEIGHT=max(Y)
			case 8
					AUTOZERO=0
					FIXEDPARAMETERS=0
					plots=1
					BIPOLAR=0
					MINWIDTH=zeros(size(peakshape))+(xx(2)-xx(1))
					delta=1
					CLIPHEIGHT=max(Y)
			case 9
					AUTOZERO=autozero
					FIXEDPARAMETERS=0
					plots=1
					BIPOLAR=0
					MINWIDTH=zeros(size(peakshape))+(xx(2)-xx(1))
					delta=1
			case 10
					AUTOZERO=autozero
					FIXEDPARAMETERS=fixedparameters
					plots=1
					BIPOLAR=0
					MINWIDTH=zeros(size(peakshape))+(xx(2)-xx(1))
					delta=1
			case 11
					AUTOZERO=autozero
					FIXEDPARAMETERS=fixedparameters
					BIPOLAR=0
					MINWIDTH=zeros(size(peakshape))+(xx(2)-xx(1))
					delta=1
					CLIPHEIGHT=max(Y)
			case 12
					AUTOZERO=autozero
					FIXEDPARAMETERS=fixedparameters
					BIPOLAR=bipolar
					MINWIDTH=zeros(size(peakshape))+(xx(2)-xx(1))
					delta=1
					CLIPHEIGHT=max(Y)
			case 13
					AUTOZERO=autozero
					FIXEDPARAMETERS=fixedparameters
					BIPOLAR=bipolar
					MINWIDTH=minwidth
					delta=1
			case 14
					AUTOZERO=autozero
					FIXEDPARAMETERS=fixedparameters
					BIPOLAR=bipolar
					MINWIDTH=minwidth
					delta=DELTA
					CLIPHEIGHT=max(Y)
			case 15
					AUTOZERO=autozero
					FIXEDPARAMETERS=fixedparameters
					BIPOLAR=bipolar
					MINWIDTH=minwidth
					delta=DELTA
					CLIPHEIGHT=clipheight
			otherwise
		# switch nargin

		# Saturation Code, skips points greater than set maximum
		if CLIPHEIGHT<max(Y),
			apnt=1
			for pnt=1:length(xx),
					if yy(pnt)<CLIPHEIGHT,
							axx(apnt)=xx(pnt)
							ayy(apnt)=yy(pnt)
							apnt=apnt+1
					
			
			xx=axxyy=ayy

		# Default values for placeholder zeros1
		if NumTrials==0NumTrials=1
		if isscalar(peakshape),
		else:
			# disp('peakshape is vector')
			shapesvector=peakshape
			NumPeaks=length(peakshape)
			peakshape=22

		if peakshape==0peakshape=1
		if NumPeaks==0NumPeaks=1
		if start==0start=calcstart(xx,NumPeaks,xoffset)
		if FIXEDPARAMETERS==0, FIXEDPARAMETERS=length(xx)/10
		if peakshape==16FIXEDPOSITIONS=fixedparameters
		if peakshape==17FIXEDPOSITIONS=fixedparameters
		if AUTOZERO>3,AUTOZERO=3,
		if AUTOZERO<0,AUTOZERO=0,
		Heights=zeros(1,NumPeaks)
		FitResults=zeros(NumPeaks,6)

		# # Remove linear baseline from data segment if AUTOZERO==1
		baseline=0
		bkgcoef=0
		bkgsize=round(length(xx)/10)
		if bkgsize<2,bkgsize=2
		lxx=length(xx)
		if AUTOZERO==1, # linear autozero operation  
			XX1=xx(1:round(lxx/bkgsize))
			XX2=xx((lxx-round(lxx/bkgsize)):lxx)
			Y1=yy(1:(round(length(xx)/bkgsize)))
			Y2=yy((lxx-round(lxx/bkgsize)):lxx)
			bkgcoef=polyfit([XX1,XX2],[Y1,Y2],1)  # Fit straight line to sub-group of points
			bkg=polyval(bkgcoef,xx)
			yy=yy-bkg
		# if
		if AUTOZERO==2, # Quadratic autozero operation  
			XX1=xx(1:round(lxx/bkgsize))
			XX2=xx((lxx-round(lxx/bkgsize)):lxx)
			Y1=yy(1:round(length(xx)/bkgsize))
			Y2=yy((lxx-round(lxx/bkgsize)):lxx)
			bkgcoef=polyfit([XX1,XX2],[Y1,Y2],2)  # Fit parabola to sub-group of points
			bkg=polyval(bkgcoef,xx)
			yy=yy-bkg
		# if autozero

		PEAKHEIGHTS=zeros(1,NumPeaks)
		n=length(xx)
		newstart=start
		# Assign ShapStrings
		switch peakshape(1)
			case 1
					ShapeString='Gaussian'
			case 2
					ShapeString='Lorentzian'
			case 3
					ShapeString='Logistic'
			case 4
					ShapeString='Pearson'
			case 5
					ShapeString='ExpGaussian'
			case 6
					ShapeString='Equal width Gaussians'
			case 7
					ShapeString='Equal width Lorentzians'
			case 8
					ShapeString='Exp. equal width Gaussians'
			case 9
					ShapeString='Exponential Pulse'
			case 10
					ShapeString='Up Sigmoid (logistic def)'
			case 23
					ShapeString='Down Sigmoid (logistic def)'  
			case 11
					ShapeString='Fixed-width Gaussian'
			case 12
					ShapeString='Fixed-width Lorentzian'
			case 13
					ShapeString='Gaussian/Lorentzian bl'
			case 14
					ShapeString='BiGaussian'    
			case 15
					ShapeString='Breit-Wigner-Fano'   
			case 16
					ShapeString='Fixed-position Gaussians'
			case 17
					ShapeString='Fixed-position Lorentzians'
			case 18
					ShapeString='Exp. Lorentzian'
			case 19
					ShapeString='Alpha def'
			case 20
					ShapeString='Voigt (equal alphas)'
			case 21
					ShapeString='triangular'
			case 22
					ShapeString=num2str(shapesvector)
			case 24
					ShapeString='Negative Binomial Distribution'
			case 25
					ShapeString='Lognormal Distribution'
			case 26
					ShapeString='slope'
			case 27
					ShapeString='First derivative'
			case 28
					ShapeString='Polynomial'
			case 29
					ShapeString='Segmented linear'
			case 30
					ShapeString='Voigt (variable alphas)'
			case 31
					ShapeString='ExpGaussian (var. time constant)'
			case 32
					ShapeString='Pearson (var. shape constant)'
			case 33
					ShapeString='Variable Gaussian/Lorentzian'
			otherwise
		# switch peakshape

		# Perform peak fitting for selected peak shape using fminsearch def
		options = optimset('TolX',.001,'Display','off','MaxFunEvals',1000 )
		LowestError=1000 # or any big number greater than largest error expected
		FitParameters=zeros(1,NumPeaks.*2) 
		BestStart=zeros(1,NumPeaks.*2) 
		height=zeros(1,NumPeaks) 
		bestmodel=zeros(size(yy))

		for k=1:NumTrials, 
			# StartMatrix(k,:)=newstart
			# disp(['Trial number ' num2str(k) ] ) # optionally prints the current trial number as progress indicator
			switch peakshape(1)
					case 1
							TrialParameters=fminsearch(@(lambda)(fitgaussian(lambda,xx,yy)),newstart,options)
							for Peak=1:NumPeaks
									if TrialParameters(2*Peak)<MINWIDTH,
											TrialParameters(2*Peak)=MINWIDTH
									
							
					case 2
							TrialParameters=fminsearch(@(lambda)(fitlorentzian(lambda,xx,yy)),newstart,options)
							for Peak=1:NumPeaks
									if TrialParameters(2*Peak)<MINWIDTH,
											TrialParameters(2*Peak)=MINWIDTH
									
							
					case 3
							TrialParameters=fminsearch(@(lambda)(fitlogistic(lambda,xx,yy)),newstart,options)
							for Peak=1:NumPeaks
									if TrialParameters(2*Peak)<MINWIDTH,
											TrialParameters(2*Peak)=MINWIDTH
									
							
					case 4
							TrialParameters=fminsearch(@(lambda)(fitpearson(lambda,xx,yy,extra)),newstart,options)
							for Peak=1:NumPeaks
									if TrialParameters(2*Peak)<MINWIDTH,
											TrialParameters(2*Peak)=MINWIDTH
									
							
					case 5
							zxx=[zeros(size(xx)) xx zeros(size(xx)) ]
							zyy=[zeros(size(yy)) yy zeros(size(yy)) ]
							TrialParameters=fminsearch(@(lambda)(fitexpgaussian(lambda,zxx,zyy,-extra)),newstart,options)
							for Peak=1:NumPeaks
									if TrialParameters(2*Peak)<MINWIDTH,
											TrialParameters(2*Peak)=MINWIDTH
									
							
					case 6
							cwnewstart(1)=newstart(1)
							for pc=2:NumPeaks,
									cwnewstart(pc)=newstart(2.*pc-1)
							
							cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5
							TrialParameters=fminsearch(@(lambda)(fitewgaussian(lambda,xx,yy)),cwnewstart,options)
							for Peak=1:NumPeaks
									if TrialParameters(NumPeaks+1)<MINWIDTH,
											TrialParameters(NumPeaks+1)=MINWIDTH
									
							
					case 7
							cwnewstart(1)=newstart(1)
							for pc=2:NumPeaks,
									cwnewstart(pc)=newstart(2.*pc-1)
							
							cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5
							TrialParameters=fminsearch(@(lambda)(fitewlorentzian(lambda,xx,yy)),cwnewstart,options)
							for Peak=1:NumPeaks
									if TrialParameters(NumPeaks+1)<MINWIDTH,
											TrialParameters(NumPeaks+1)=MINWIDTH
									
							
					case 8
							cwnewstart(1)=newstart(1)
							for pc=2:NumPeaks,
									cwnewstart(pc)=newstart(2.*pc-1)
							
							cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5
							TrialParameters=fminsearch(@(lambda)(fitexpewgaussian(lambda,xx,yy,-extra)),cwnewstart,options)
							for Peak=1:NumPeaks
									if TrialParameters(NumPeaks+1)<MINWIDTH,
											TrialParameters(NumPeaks+1)=MINWIDTH
									
							
					case 9
							TrialParameters=fminsearch(@(lambda)(fitexppulse(lambda,xx,yy)),newstart,options)
							for Peak=1:NumPeaks
									if TrialParameters(2*Peak)<MINWIDTH,
											TrialParameters(2*Peak)=MINWIDTH
									
							
					case 10
							TrialParameters=fminsearch(@(lambda)(fitupsigmoid(lambda,xx,yy)),newstart,options)
							for Peak=1:NumPeaks
									if TrialParameters(2*Peak)<MINWIDTH,
											TrialParameters(2*Peak)=MINWIDTH
									
							
					case 23
							TrialParameters=fminsearch(@(lambda)(fitdownsigmoid(lambda,xx,yy)),newstart,options)
							for Peak=1:NumPeaks
									if TrialParameters(2*Peak)<MINWIDTH,
											TrialParameters(2*Peak)=MINWIDTH
									
							
					case 11
							fixedstart=[]
							for pc=1:NumPeaks,
									fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1)
							
							TrialParameters=fminsearch(@(lambda)(FitFWGaussian(lambda,xx,yy)),fixedstart,options)
					case 12
							fixedstart=[]
							for pc=1:NumPeaks,
									fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1)
							
							TrialParameters=fminsearch(@(lambda)(FitFWLorentzian(lambda,xx,yy)),fixedstart,options)
					case 13
							TrialParameters=fminsearch(@(lambda)(fitGL(lambda,xx,yy,extra)),newstart,options)
							for Peak=1:NumPeaks
									if TrialParameters(2*Peak)<MINWIDTH,
											TrialParameters(2*Peak)=MINWIDTH
									
							
					case 14
							TrialParameters=fminsearch(@(lambda)(fitBiGaussian(lambda,xx,yy,extra)),newstart,options)
							for Peak=1:NumPeaks
									if TrialParameters(2*Peak)<MINWIDTH,
											TrialParameters(2*Peak)=MINWIDTH
									
							
					case 15
							TrialParameters=fminsearch(@(lambda)(fitBWF(lambda,xx,yy,extra)),newstart,options)
							for Peak=1:NumPeaks
									if TrialParameters(2*Peak)<MINWIDTH,
											TrialParameters(2*Peak)=MINWIDTH
									
							
					case 16
							fixedstart=[]
							for pc=1:NumPeaks,
									fixedstart(pc)=(max(xx)-min(xx))./(NumPeaks+1)
									fixedstart(pc)=fixedstart(pc)+.1*(rand-.5).*fixedstart(pc)
							
							TrialParameters=fminsearch(@(lambda)(FitFPGaussian(lambda,xx,yy)),fixedstart,options)
							for Peak=1:NumPeaks
									if TrialParameters(Peak)<MINWIDTH,
											TrialParameters(Peak)=MINWIDTH
									
							
					case 17
							fixedstart=[]
							for pc=1:NumPeaks,
									fixedstart(pc)=(max(xx)-min(xx))./(NumPeaks+1)
									fixedstart(pc)=fixedstart(pc)+.1*(rand-.5).*fixedstart(pc)
							
							TrialParameters=fminsearch(@(lambda)(FitFPLorentzian(lambda,xx,yy)),fixedstart,options)
							for Peak=1:NumPeaks
									if TrialParameters(Peak)<MINWIDTH,
											TrialParameters(Peak)=MINWIDTH
									
							
					case 18
							zxx=[zeros(size(xx)) xx zeros(size(xx)) ]
							zyy=[ones(size(yy)).*yy(1) yy zeros(size(yy)).*yy(length(yy)) ]
							TrialParameters=fminsearch(@(lambda)(fitexplorentzian(lambda,zxx,zyy,-extra)),newstart,options)
							for Peak=1:NumPeaks
									if TrialParameters(2*Peak)<MINWIDTH,
											TrialParameters(2*Peak)=MINWIDTH
									
							
					case 19
							TrialParameters=fminsearch(@(lambda)(fitalphadef(lambda,xx,yy)),newstart,options)
							for Peak=1:NumPeaks
									if TrialParameters(2*Peak)<MINWIDTH,
											TrialParameters(2*Peak)=MINWIDTH
									
							
					case 20
							TrialParameters=fminsearch(@(lambda)(fitvoigt(lambda,xx,yy,extra)),newstart,options)
							for Peak=1:NumPeaks
									if TrialParameters(2*Peak)<MINWIDTH,
											TrialParameters(2*Peak)=MINWIDTH
									
							
					case 21
							TrialParameters=fminsearch(@(lambda)(fittriangular(lambda,xx,yy)),newstart,options)
							for Peak=1:NumPeaks
									if TrialParameters(2*Peak)<MINWIDTH,
											TrialParameters(2*Peak)=MINWIDTH
									
							
					case 22
							TrialParameters=fminsearch(@(lambda)(fitmultiple(lambda,xx,yy,shapesvector,extra)),newstart,options)
							for Peak=1:NumPeaks
									if TrialParameters(2*Peak)<MINWIDTH(Peak),
											TrialParameters(2*Peak)=MINWIDTH(Peak)
									
							
					case 24
							TrialParameters=fminsearch(@(lambda)(fitnbinpdf(lambda,xx,yy)),newstart,options)
							for Peak=1:NumPeaks
									if TrialParameters(2*Peak)<MINWIDTH,
											TrialParameters(2*Peak)=MINWIDTH
									
							
					case 25
							TrialParameters=fminsearch(@(lambda)(fitlognpdf(lambda,xx,yy)),newstart,options)
							for Peak=1:NumPeaks
									if TrialParameters(2*Peak)<MINWIDTH,
											TrialParameters(2*Peak)=MINWIDTH
									
							
					case 26
							TrialParameters=fminsearch(@(lambda)(fitlinslope(lambda,xx,yy)),polyfit(xx,yy,1),options)
							 coeff=TrialParameters
					case 27
							TrialParameters=fminsearch(@(lambda)(fitd1gauss(lambda,xx,yy)),newstart,options)
					case 28
							coeff=fitpolynomial(xx,yy,extra)
							TrialParameters=coeff
					case 29
							cnewstart(1)=newstart(1)
							for pc=2:NumPeaks,
									cnewstart(pc)=newstart(2.*pc-1)+(delta*(rand-.5)/50)
							
							TrialParameters=fminsearch(@(lambda)(fitsegmented(lambda,xx,yy)),cnewstart,options)
					case 30
							nn=max(xx)-min(xx)
							start=[]
							startpos=[nn/(NumPeaks+1):nn/(NumPeaks+1):nn-(nn/(NumPeaks+1))]+min(xx)
							for marker=1:NumPeaks,
									markx=startpos(marker)+ xoffset
									start=[start markx nn/5 extra]
							 # for marker
							 newstart=start
							for parameter=1:3:3*NumPeaks,
									newstart(parameter)=newstart(parameter)*(1+randn/100)
									newstart(parameter+1)=newstart(parameter+1)*(1+randn/20)
									newstart(parameter+2)=newstart(parameter+1)*(1+randn/20)
							
							TrialParameters=fminsearch(@(lambda)(fitvoigtv(lambda,xx,yy)),newstart)
					 case 31
							nn=max(xx)-min(xx)
							start=[]
							startpos=[nn/(NumPeaks+1):nn/(NumPeaks+1):nn-(nn/(NumPeaks+1))]+min(xx)
							for marker=1:NumPeaks,
									markx=startpos(marker)+ xoffset
									start=[start markx nn/5 extra]
							 # for marker
							 newstart=start
							for parameter=1:3:3*NumPeaks,
									newstart(parameter)=newstart(parameter)*(1+randn/100)
									newstart(parameter+1)=newstart(parameter+1)*(1+randn/20)
									newstart(parameter+2)=newstart(parameter+1)*(1+randn/20)
							
							# newstart=newstart
							zxx=[zeros(size(xx)) xx zeros(size(xx)) ]
							zyy=[ones(size(yy)).*yy(1) yy zeros(size(yy)).*yy(length(yy)) ]
							TrialParameters=fminsearch(@(lambda)(fitexpgaussianv(lambda,zxx,zyy)),newstart)
					case 32
							nn=max(xx)-min(xx)
							start=[]
							startpos=[nn/(NumPeaks+1):nn/(NumPeaks+1):nn-(nn/(NumPeaks+1))]+min(xx)
							for marker=1:NumPeaks,
									markx=startpos(marker)+ xoffset
									start=[start markx nn/5 extra]
							 # for marker
							 newstart=start
							for parameter=1:3:3*NumPeaks,
									newstart(parameter)=newstart(parameter)*(1+randn/100)
									newstart(parameter+1)=newstart(parameter+1)*(1+randn/20)
									newstart(parameter+2)=newstart(parameter+1)*(1+randn/20)
							
							# newstart=newstart
							TrialParameters=fminsearch(@(lambda)(fitpearsonv(lambda,xx,yy)),newstart)
					case 33
							 nn=max(xx)-min(xx)
							start=[]
							startpos=[nn/(NumPeaks+1):nn/(NumPeaks+1):nn-(nn/(NumPeaks+1))]+min(xx)
							for marker=1:NumPeaks,
									markx=startpos(marker)+ xoffset
									start=[start markx nn/5 extra]
							 # for marker
							 newstart=start
							for parameter=1:3:3*NumPeaks,
									newstart(parameter)=newstart(parameter)*(1+randn/100)
									newstart(parameter+1)=newstart(parameter+1)*(1+randn/20)
									newstart(parameter+2)=newstart(parameter+1)*(1+randn/20)
							
							# newstart=newstart
							TrialParameters=fminsearch(@(lambda)(fitGLv(lambda,xx,yy)),newstart)
					otherwise
			 # switch peakshape

		# Construct model from Trial parameters
		A=zeros(NumPeaks,n)
		for m=1:NumPeaks,
			switch peakshape(1)
					case 1
							A(m,:)=gaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m))
					case 2
							A(m,:)=lorentzian(xx,TrialParameters(2*m-1),TrialParameters(2*m))
					case 3
							A(m,:)=logistic(xx,TrialParameters(2*m-1),TrialParameters(2*m))
					case 4
							A(m,:)=pearson(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra)
					case 5
							A(m,:)=expgaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m),-extra)'
					case 6
							A(m,:)=gaussian(xx,TrialParameters(m),TrialParameters(NumPeaks+1))
					case 7
							A(m,:)=lorentzian(xx,TrialParameters(m),TrialParameters(NumPeaks+1))
					case 8
							A(m,:)=expgaussian(xx,TrialParameters(m),TrialParameters(NumPeaks+1),-extra)'
					case 9
							A(m,:)=exppulse(xx,TrialParameters(2*m-1),TrialParameters(2*m))
					case 10
							A(m,:)=upsigmoid(xx,TrialParameters(2*m-1),TrialParameters(2*m))
					case 11
							A(m,:)=gaussian(xx,TrialParameters(m),FIXEDPARAMETERS)
					case 12
							A(m,:)=lorentzian(xx,TrialParameters(m),FIXEDPARAMETERS)
					case 13
							A(m,:)=GL(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra)
					case 14
							A(m,:)=BiGaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra)
					case 15
							A(m,:)=BWF(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra)        
					case 16
							A(m,:)=gaussian(xx,FIXEDPOSITIONS(m),TrialParameters(m))
					case 17
							A(m,:)=lorentzian(xx,FIXEDPOSITIONS(m),TrialParameters(m))
					case 18
							A(m,:)=explorentzian(xx,TrialParameters(2*m-1),TrialParameters(2*m),-extra)'
					case 19
							A(m,:)=alphadef(xx,TrialParameters(2*m-1),TrialParameters(2*m))
					case 20
							A(m,:)=voigt(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra)        
					case 21
							A(m,:)=triangular(xx,TrialParameters(2*m-1),TrialParameters(2*m))
					case 22
							A(m,:)=peakdef(shapesvector(m),xx,TrialParameters(2*m-1),TrialParameters(2*m),extra(m))        
					case 23
							A(m,:)=downsigmoid(xx,TrialParameters(2*m-1),TrialParameters(2*m))        
					case 24
							A(m,:)=nbinpdf(xx,TrialParameters(2*m-1),TrialParameters(2*m))
					case 25
							A(m,:)=lognormal(xx,TrialParameters(2*m-1),TrialParameters(2*m))
					case 26
							A(m,:)=linslope(xx,TrialParameters(2*m-1),TrialParameters(2*m))
					case 27
							A(m,:)=d1gauss(xx,TrialParameters(2*m-1),TrialParameters(2*m))
					case 28
							A(m,:)=polynomial(xx,coeff)
					case 29
							A(m,:)=segmented(xx,yy,PEAKHEIGHTS)
					case 30
							A(m,:)=voigt(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m))        
					case 31
							A(m,:)=expgaussian(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),-TrialParameters(3*m))        
					case 32
							A(m,:)=pearson(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m))        
					case 33
							A(m,:)=GL(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m))        
					otherwise
			 # switch
			for parameter=1:2:2*NumPeaks,
					newstart(parameter)=newstart(parameter)*(1+delta*(rand-.5)/500)
					newstart(parameter+1)=newstart(parameter+1)*(1+delta*(rand-.5)/100)
			
		# for NumPeaks

		# Multiplies each row by the corresponding amplitude and adds them up
		if peakshape(1)==29, # Segmented linear
			model=segmented(xx,yy,PEAKHEIGHTS)
			TrialParameters=PEAKHEIGHTS
			Heights=ones(size(PEAKHEIGHTS))
		else:
			if AUTOZERO==3,
					baseline=PEAKHEIGHTS(1)
					Heights=PEAKHEIGHTS(2:1+NumPeaks)
					model=Heights'*A+baseline
			else:
		#          size(PEAKHEIGHTS) # error check
		#          size(A)
					model=PEAKHEIGHTS'*A
					Heights=PEAKHEIGHTS
					baseline=0
			

		if peakshape(1)==28, # polynomial
			model=polynomial(xx,coeff)
			TrialParameters=PEAKHEIGHTS
			Heights=ones(size(PEAKHEIGHTS))

		# Compare trial model to data segment and compute the fit error
			MeanFitError=100*norm(yy-model)./(sqrt(n)*max(yy))
		# Take only the single fit that has the lowest MeanFitError
		if MeanFitError<LowestError, 
				if min(Heights)>=-BIPOLAR*10^100,  # Consider only fits with positive peak heights
					LowestError=MeanFitError  # Assign LowestError to the lowest MeanFitError
					FitParameters=TrialParameters  # Assign FitParameters to the fit with the lowest MeanFitError
					BestStart=newstart # Assign BestStart to the start with the lowest MeanFitError
					height=Heights # Assign height to the PEAKHEIGHTS with the lowest MeanFitError
					bestmodel=model # Assign bestmodel to the model with the lowest MeanFitError
				 # if min(PEAKHEIGHTS)>0
		 # if MeanFitError<LowestError
		#  ErrorVector(k)=MeanFitError
		# for k (NumTrials)
			Rsquared=1-(norm(yy-bestmodel)./norm(yy-mean(yy)))
			SStot=sum((yy-mean(yy)).^2)
			SSres=sum((yy-bestmodel).^2)
			Rsquared=1-(SSres./SStot)
			GOF=[LowestError Rsquared]
		# Uncomment following 4 lines to monitor trail fit starts and errors.
		# StartMatrix=StartMatrix
		# ErrorVector=ErrorVector
		# matrix=[StartMatrix ErrorVector']
		# std(StartMatrix)
		# Construct model from best-fit parameters
		AA=zeros(NumPeaks,600)
		xxx=linspace(min(xx),max(xx),600)
		# xxx=linspace(min(xx)-length(xx),max(xx)+length(xx),200)
		for m=1:NumPeaks,
		 switch peakshape(1)
			case 1
					AA(m,:)=gaussian(xxx,FitParameters(2*m-1),FitParameters(2*m))
			case 2
					AA(m,:)=lorentzian(xxx,FitParameters(2*m-1),FitParameters(2*m))
			case 3
					AA(m,:)=logistic(xxx,FitParameters(2*m-1),FitParameters(2*m))
			case 4
					AA(m,:)=pearson(xxx,FitParameters(2*m-1),FitParameters(2*m),extra)
			case 5
					AA(m,:)=expgaussian(xxx,FitParameters(2*m-1),FitParameters(2*m),-extra*length(xxx)./length(xx))'
			case 6
					AA(m,:)=gaussian(xxx,FitParameters(m),FitParameters(NumPeaks+1))
			case 7
					AA(m,:)=lorentzian(xxx,FitParameters(m),FitParameters(NumPeaks+1))
			case 8
					AA(m,:)=expgaussian(xxx,FitParameters(m),FitParameters(NumPeaks+1),-extra*length(xxx)./length(xx))'
			case 9
					AA(m,:)=exppulse(xxx,FitParameters(2*m-1),FitParameters(2*m))  
			case 10
					AA(m,:)=upsigmoid(xxx,FitParameters(2*m-1),FitParameters(2*m))   
			case 11
					AA(m,:)=gaussian(xxx,FitParameters(m),FIXEDPARAMETERS)
			case 12
					AA(m,:)=lorentzian(xxx,FitParameters(m),FIXEDPARAMETERS)
			case 13
					AA(m,:)=GL(xxx,FitParameters(2*m-1),FitParameters(2*m),extra)
			case 14
					AA(m,:)=BiGaussian(xxx,FitParameters(2*m-1),FitParameters(2*m),extra)       
			case 15
					AA(m,:)=BWF(xxx,FitParameters(2*m-1),FitParameters(2*m),extra)       
			case 16
					AA(m,:)=gaussian(xxx,FIXEDPOSITIONS(m),FitParameters(m))
			case 17
					AA(m,:)=lorentzian(xxx,FIXEDPOSITIONS(m),FitParameters(m))
			case 18
					AA(m,:)=explorentzian(xxx,FitParameters(2*m-1),FitParameters(2*m),-extra*length(xxx)./length(xx))'
			case 19
					AA(m,:)=alphadef(xxx,FitParameters(2*m-1),FitParameters(2*m))
			case 20
					AA(m,:)=voigt(xxx,FitParameters(2*m-1),FitParameters(2*m),extra)       
			case 21
					AA(m,:)=triangular(xxx,FitParameters(2*m-1),FitParameters(2*m))
			case 22
					AA(m,:)=peakdef(shapesvector(m),xxx,FitParameters(2*m-1),FitParameters(2*m),extra(m))        
			case 23
					AA(m,:)=downsigmoid(xxx,FitParameters(2*m-1),FitParameters(2*m))  
			case 24
					AA(m,:)=nbinpdf(xxx,FitParameters(2*m-1),FitParameters(2*m))    
			case 25
					AA(m,:)=lognormal(xxx,FitParameters(2*m-1),FitParameters(2*m))    
			case 26
					AA(m,:)=linslope(xxx,FitParameters(2*m-1),FitParameters(2*m))   
			case 27
					AA(m,:)=d1gauss(xxx,FitParameters(2*m-1),FitParameters(2*m))  
			case 28
					AA(m,:)=polynomial(xxx,coeff)
			case 29
			case 30
					AA(m,:)=voigt(xxx,FitParameters(3*m-2),FitParameters(3*m-1),FitParameters(3*m))        
			case 31
					AA(m,:)=expgaussian(xxx,FitParameters(3*m-2),FitParameters(3*m-1),-FitParameters(3*m)*length(xxx)./length(xx))        
			case 32
					AA(m,:)=pearson(xxx,FitParameters(3*m-2),FitParameters(3*m-1),FitParameters(3*m))        
			case 33
					AA(m,:)=GL(xxx,FitParameters(3*m-2),FitParameters(3*m-1),FitParameters(3*m)) 
				 otherwise
		 # switch
		# for NumPeaks

		# Multiplies each row by the corresponding amplitude and adds them up
		if peakshape(1)==29, # Segmented linear
			mmodel=segmented(xx,yy,PEAKHEIGHTS)
			baseline=0
		else:
			heightsize=size(height')
			AAsize=size(AA)
			if heightsize(2)==AAsize(1),
					mmodel=height'*AA+baseline
			else:
					mmodel=height*AA+baseline
			

		# Top half of the figure shows original signal and the fitted model.
		if plots,
			subplot(2,1,1)plot(xx+xoffset,yy,'b.') # Plot the original signal in blue dots
			hold on

		if peakshape(1)==28, # Polynomial
			 yi=polynomial(xxx,coeff)
		else:
			for m=1:NumPeaks,
					if plots, plot(xxx+xoffset,height(m)*AA(m,:)+baseline,'g'),  # Plot the individual component peaks in green lines
					area(m)=trapz(xxx+xoffset,height(m)*AA(m,:)) # Compute the area of each component peak using trapezoidal method
					yi(m,:)=height(m)*AA(m,:) # Place y values of individual model peaks into matrix yi
			

		xi=xxx+xoffset # Place the x-values of the individual model peaks into xi

		if plots,
			# Mark starting peak positions with vertical dashed magenta lines
			if peakshape(1)==16||peakshape(1)==17
			else:
					if peakshape(1)==29, # Segmented linear
							subplot(2,1,1)plot([PEAKHEIGHTS' PEAKHEIGHTS'],[0 max(yy)],'m--')
					else:
							for marker=1:NumPeaks,
									markx=BestStart((2*marker)-1)
									subplot(2,1,1)plot([markx+xoffset markx+xoffset],[0 max(yy)],'m--')
							 # for
					
			 # if peakshape

			# Plot the total model (sum of component peaks) in red lines
			if peakshape(1)==29, # Segmented linear
					mmodel=segmented(xx,yy,PEAKHEIGHTS)
				 plot(xx+xoffset,mmodel,'r')  
			else:
				 plot(xxx+xoffset,mmodel,'r')  
			
			hold off
			lyy=min(yy)
			uyy=max(yy)+(max(yy)-min(yy))/10
			if BIPOLAR,
					axis([min(xx) max(xx) lyy uyy])
					ylabel('+ - mode')
			else:
					axis([min(xx) max(xx) 0 uyy])
					ylabel('+ mode')
			
			switch AUTOZERO,
					case 0
							title(['peakfit.m Version 7   No baseline correction'])
					case 1
							title(['peakfit.m Version 7   Linear baseline subtraction'])
					case 2
							title(['peakfit.m Version 7   Quadratic subtraction baseline'])
					case 3
							title(['peakfit.m Version 7   Flat baseline correction'])
			

			switch peakshape(1)
					case {4,20}
							xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Min. Width = ' num2str(MINWIDTH) '      Shape Constant = ' num2str(extra)  '     Error = ' num2str(round(1000*LowestError)/1000) '#   R2 = ' num2str(round(100000*Rsquared)/100000) ] )
					case {5,8,18}
							xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Min. Width = ' num2str(MINWIDTH) '      Time Constant = ' num2str(extra)   '     Error = ' num2str(round(1000*LowestError)/1000) '#   R2 = ' num2str(round(100000*Rsquared)/100000)  ] )
					case 13
							xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Min. Width = ' num2str(MINWIDTH) '      # Gaussian = ' num2str(extra)   '     Error = ' num2str(round(1000*LowestError)/1000) '#   R2 = ' num2str(round(100000*Rsquared)/100000)  ] )
					case {14,15,22}
							xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Min. Width = ' num2str(MINWIDTH) '      extra = ' num2str(extra)  '     Error = ' num2str(round(1000*LowestError)/1000) '#   R2 = ' num2str(round(100000*Rsquared)/100000) ] )
					case 28
							xlabel(['Shape = ' ShapeString '      Order = ' num2str(extra)  '     Error = ' num2str(round(1000*LowestError)/1000) '#  R2 = ' num2str(round(1000*LowestError)/1000) ] )
					otherwise
							if peakshape(1)==29, # Segmented linear
									xlabel(['Breakpoints = ' num2str(NumPeaks) '     Shape = ' ShapeString  '     Error = ' num2str(round(1000*LowestError)/1000) '#  R2 = ' num2str(round(100000*Rsquared)/100000) ] )
							else:
									xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Min. Width = ' num2str(MINWIDTH)  '     Error = ' num2str(round(1000*LowestError)/1000) '#  R2 = ' num2str(round(100000*Rsquared)/100000) ] )
							 # if peakshape(1)==29
			 # switch peakshape(1)

			# Bottom half of the figure shows the residuals and displays RMS error
			# between original signal and model
			residual=yy-bestmodel
			subplot(2,1,2)plot(xx+xoffset,residual,'r.')
			axis([min(xx)+xoffset max(xx)+xoffset min(residual) max(residual)])
			xlabel('Residual Plot')
			if NumTrials>1,
				 title(['Best of ' num2str(NumTrials) ' fits'])
			else:
				 title(['Single fit'])
			
		# if plots

		# Put results into a matrix FitResults, one row for each peak, showing peak index number,
		# position, amplitude, and width.
		FitResults=zeros(NumPeaks,6)
		switch peakshape(1),
			case {6,7,8}, # equal-width peak models only
					for m=1:NumPeaks,
							if m==1,
									FitResults=[[round(m) FitParameters(m)+xoffset height(m) abs(FitParameters(NumPeaks+1)) area(m)]]
							else:
									FitResults=[FitResults  [round(m) FitParameters(m)+xoffset height(m) abs(FitParameters(NumPeaks+1)) area(m)]]
							
					
			case {11,12}, # Fixed-width shapes only
					for m=1:NumPeaks,
							if m==1,
									FitResults=[[round(m) FitParameters(m)+xoffset height(m) FIXEDPARAMETERS area(m)]]
							else:
									FitResults=[FitResults  [round(m) FitParameters(m)+xoffset height(m) FIXEDPARAMETERS area(m)]]
							
					
			case {16,17}, # Fixed-position shapes only
					for m=1:NumPeaks,
							if m==1,
									FitResults=[round(m) FIXEDPOSITIONS(m) height(m) FitParameters(m) area(m)]
							else:
									FitResults=[FitResults  [round(m) FIXEDPOSITIONS(m) height(m) FitParameters(m) area(m)]]
							
					
			case 28,   # Simple polynomial fit
					FitResults=PEAKHEIGHTS
			case 29, # Segmented linear fit
					FitResults=PEAKHEIGHTS
			case {30,31,32,33} # Special case of shapes with 3 iterated variables
					for m=1:NumPeaks,
							if m==1,
									FitResults=[round(m) FitParameters(3*m-2) height(m) abs(FitParameters(3*m-1)) area(m) FitParameters(3*m)]
							else:
									FitResults=[FitResults  [round(m) FitParameters(3*m-2) height(m) abs(FitParameters(3*m-1)) area(m)] FitParameters(3*m)]
							
					
			otherwise # Normal shapes with 2 iterated variables
					for m=1:NumPeaks,
							if m==1,
									FitResults=[round(m) FitParameters(2*m-1)+xoffset height(m) abs(FitParameters(2*m)) area(m)]
							else:
									FitResults=[FitResults  [round(m) FitParameters(2*m-1)+xoffset height(m) abs(FitParameters(2*m)) area(m)]]
							 # if m==1
					 # for m=1:NumPeaks,
		# switch peakshape(1)

		# Display Fit Results on lower graph
		if plots,
			# Display Fit Results on lower  graph
			subplot(2,1,2)
			startx=min(xx)+(max(xx)-min(xx))./20
			dxx=(max(xx)-min(xx))./10
			dyy=((max(residual)-min(residual))./10)
			starty=max(residual)-dyy
			FigureSize=get(gcf,'Position')
			switch peakshape(1)
					case {9,19,11}  # Pulse and sigmoid shapes only
							text(startx,starty+dyy/2,['Peak #          tau1           Height           tau2             Area'] )
					case 28, # Polynomial
							text(startx,starty+dyy/2,['Polynomial coefficients'] )
					case 29 # Segmented linear
							 text(startx,starty+dyy/2,['x-axis breakpoints'] )
					case {30,31,32,33} # Special case of shapes with 3 iterated variables
							text(startx,starty+dyy/2,['Peak #          Position        Height         Width             Area       Shape factor'] )            
					otherwise
							text(startx,starty+dyy/2,['Peak #          Position        Height         Width             Area '] )
			
			# Display FitResults using sprintf
			if peakshape(1)==28||peakshape(1)==29, # Polynomial or segmented linear
					for number=1:length(FitResults),
							column=1
							itemstring=sprintf('#0.4g',FitResults(number))
							xposition=startx+(1.7.*dxx.*(column-1).*(600./FigureSize(3)))
							yposition=starty-number.*dyy.*(400./FigureSize(4))
							text(xposition,yposition,['                ' itemstring])
					
			else:
					for peaknumber=1:NumPeaks,
							for column=1:5,
									itemstring=sprintf('#0.4g',FitResults(peaknumber,column))
									xposition=startx+(1.7.*dxx.*(column-1).*(600./FigureSize(3)))
									yposition=starty-peaknumber.*dyy.*(400./FigureSize(4))
									text(xposition,yposition,itemstring)
							
					
					xposition=startx
					yposition=starty-(peaknumber+1).*dyy.*(400./FigureSize(4))
					if AUTOZERO==3,
							text(xposition,yposition,[ 'Baseline= ' num2str(baseline) ])
					 # if AUTOZERO
			 # if peakshape(1)
			if peakshape(1)==30 || peakshape(1)==31 || peakshape(1)==32 || peakshape(1)==33,
					for peaknumber=1:NumPeaks,
							column=6
							itemstring=sprintf('#0.4g',FitParameters(3*peaknumber))
							xposition=startx+(1.7.*dxx.*(column-1).*(600./FigureSize(3)))
							yposition=starty-peaknumber.*dyy.*(400./FigureSize(4))
							text(xposition,yposition,itemstring)
					
			
		# if plots

		if NumArgOut==8,
			if plots,disp('Computing bootstrap sampling statistics.....'),
			BootstrapResultsMatrix=zeros(6,100,NumPeaks)
			BootstrapErrorMatrix=zeros(1,100,NumPeaks)
			clear bx by
			tic
			for trial=1:100,
					n=1
					bx=xx
					by=yy
					while n<length(xx)-1,
							if rand>.5,
									bx(n)=xx(n+1)
									by(n)=yy(n+1)
							
							n=n+1
					
					bx=bx+xoffset
					[FitResults,BootFitError]=fitpeaks(bx,by,NumPeaks,peakshape,extra,NumTrials,start,AUTOZERO,FIXEDPARAMETERS)
					for peak=1:NumPeaks,
							switch peakshape(1)
									case {30,31,32,33}
											BootstrapResultsMatrix(1:6,trial,peak)=FitResults(peak,1:6)
									otherwise
											BootstrapResultsMatrix(1:5,trial,peak)=FitResults(peak,1:5)
							
							BootstrapErrorMatrix(:,trial,peak)=BootFitError
					
			
			if plots,toc
			for peak=1:NumPeaks,
					if plots,
							disp(' ')
							disp(['Peak #',num2str(peak) '         Position    Height       Width       Area      Shape Factor'])
					 # if plots
					BootstrapMean=mean(real(BootstrapResultsMatrix(:,:,peak)'))
					BootstrapSTD=std(BootstrapResultsMatrix(:,:,peak)')
					BootstrapIQR=iqr(BootstrapResultsMatrix(:,:,peak)')
					PercentRSD=100.*BootstrapSTD./BootstrapMean
					PercentIQR=100.*BootstrapIQR./BootstrapMean
					BootstrapMean=BootstrapMean(2:6)
					BootstrapSTD=BootstrapSTD(2:6)
					BootstrapIQR=BootstrapIQR(2:6)
					PercentRSD=PercentRSD(2:6)
					PercentIQR=PercentIQR(2:6)
					if plots,
							disp(['Bootstrap Mean: ', num2str(BootstrapMean)])
							disp(['Bootstrap STD:  ', num2str(BootstrapSTD)])
							disp(['Bootstrap IQR:  ', num2str(BootstrapIQR)])
							disp(['Percent RSD:    ', num2str(PercentRSD)])
							disp(['Percent IQR:    ', num2str(PercentIQR)])
					 # if plots
					BootResults(peak,:)=[BootstrapMean BootstrapSTD PercentRSD BootstrapIQR PercentIQR]
			 # peak=1:NumPeaks,
		# if NumArgOut==8,
		if AUTOZERO==3
		else:
			baseline=bkgcoef

		# ----------------------------------------------------------------------
		def [FitResults,LowestError]=fitpeaks(xx,yy,NumPeaks,peakshape,extra,NumTrials,start,AUTOZERO,fixedparameters)
		# Based on peakfit Version 3: June, 2012. 
		global PEAKHEIGHTS FIXEDPARAMETERS AUTOZERO BIPOLAR MINWIDTH coeff
		format short g
		format compact
		warning off all
		FIXEDPARAMETERS=fixedparameters
		xoffset=0
		if start==0start=calcstart(xx,NumPeaks,xoffset)
		PEAKHEIGHTS=zeros(1,NumPeaks)
		n=length(xx)
		newstart=start
		coeff=0
		LOGPLOT=0

		# Perform peak fitting for selected peak shape using fminsearch def
		options = optimset('TolX',.001,'Display','off','MaxFunEvals',1000 )
		LowestError=1000 # or any big number greater than largest error expected
		FitParameters=zeros(1,NumPeaks.*2) 
		BestStart=zeros(1,NumPeaks.*2) 
		height=zeros(1,NumPeaks) 
		bestmodel=zeros(size(yy))
		for k=1:NumTrials,
			# StartVector=newstart
			switch peakshape(1)
					case 1
							TrialParameters=fminsearch(@(lambda)(fitgaussian(lambda,xx,yy)),newstart,options)
							for Peak=1:NumPeaks
									if TrialParameters(2*Peak)<MINWIDTH,
											TrialParameters(2*Peak)=MINWIDTH
									
							
					case 2
							TrialParameters=fminsearch(@(lambda)(fitlorentzian(lambda,xx,yy)),newstart,options)
							for Peak=1:NumPeaks
									if TrialParameters(2*Peak)<MINWIDTH,
											TrialParameters(2*Peak)=MINWIDTH
									
							
					case 3
							TrialParameters=fminsearch(@(lambda)(fitlogistic(lambda,xx,yy)),newstart,options)
							for Peak=1:NumPeaks
									if TrialParameters(2*Peak)<MINWIDTH,
											TrialParameters(2*Peak)=MINWIDTH
									
							
					case 4
							TrialParameters=fminsearch(@(lambda)(fitpearson(lambda,xx,yy,extra)),newstart,options)
							for Peak=1:NumPeaks
									if TrialParameters(2*Peak)<MINWIDTH,
											TrialParameters(2*Peak)=MINWIDTH
									
							
					case 5
							zxx=[zeros(size(xx)) xx zeros(size(xx)) ]
							zyy=[zeros(size(yy)) yy zeros(size(yy)) ]
							TrialParameters=fminsearch(@(lambda)(fitexpgaussian(lambda,zxx,zyy,-extra)),newstart,options)
							for Peak=1:NumPeaks
									if TrialParameters(2*Peak)<MINWIDTH,
											TrialParameters(2*Peak)=MINWIDTH
									
							
					case 6
							cwnewstart(1)=newstart(1)
							for pc=2:NumPeaks,
									cwnewstart(pc)=newstart(2.*pc-1)
							
							cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5
							TrialParameters=fminsearch(@(lambda)(fitewgaussian(lambda,xx,yy)),cwnewstart,options)
							for Peak=1:NumPeaks
									if TrialParameters(NumPeaks+1)<MINWIDTH,
											TrialParameters(NumPeaks+1)=MINWIDTH
									
							
					case 7
							cwnewstart(1)=newstart(1)
							for pc=2:NumPeaks,
									cwnewstart(pc)=newstart(2.*pc-1)
							
							cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5
							TrialParameters=fminsearch(@(lambda)(fitewlorentzian(lambda,xx,yy)),cwnewstart,options)
							for Peak=1:NumPeaks
									if TrialParameters(NumPeaks+1)<MINWIDTH,
											TrialParameters(NumPeaks+1)=MINWIDTH
									
							
					case 8
							cwnewstart(1)=newstart(1)
							for pc=2:NumPeaks,
									cwnewstart(pc)=newstart(2.*pc-1)
							
							cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5
							TrialParameters=fminsearch(@(lambda)(fitexpewgaussian(lambda,xx,yy,-extra)),cwnewstart,options)
							for Peak=1:NumPeaks
									if TrialParameters(NumPeaks+1)<MINWIDTH,
											TrialParameters(NumPeaks+1)=MINWIDTH
									
							
					case 9
							TrialParameters=fminsearch(@(lambda)(fitexppulse(lambda,xx,yy)),newstart,options)
							for Peak=1:NumPeaks
									if TrialParameters(2*Peak)<MINWIDTH,
											TrialParameters(2*Peak)=MINWIDTH
									
							
					case 10
							TrialParameters=fminsearch(@(lambda)(fitupsigmoid(lambda,xx,yy)),newstar,optionst)
							for Peak=1:NumPeaks
									if TrialParameters(2*Peak)<MINWIDTH,
											TrialParameters(2*Peak)=MINWIDTH
									
							
					case 11
							fixedstart=[]
							for pc=1:NumPeaks,
									fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1)
							
							TrialParameters=fminsearch(@(lambda)(FitFWGaussian(lambda,xx,yy)),fixedstart,options)
							for Peak=1:NumPeaks
									if TrialParameters(2*Peak)<MINWIDTH,
											TrialParameters(2*Peak)=MINWIDTH
									
							
					case 12
							fixedstart=[]
							for pc=1:NumPeaks,
									fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1)
							
							TrialParameters=fminsearch(@(lambda)(FitFWLorentzian(lambda,xx,yy)),fixedstart,options)
							for Peak=1:NumPeaks
									if TrialParameters(2*Peak)<MINWIDTH,
											TrialParameters(2*Peak)=MINWIDTH
									
							
					case 13
							TrialParameters=fminsearch(@(lambda)(fitGL(lambda,xx,yy,extra)),newstart,options)
							for Peak=1:NumPeaks
									if TrialParameters(2*Peak)<MINWIDTH,
											TrialParameters(2*Peak)=MINWIDTH
									
							
					case 14
							TrialParameters=fminsearch(@(lambda)(fitBiGaussian(lambda,xx,yy,extra)),newstart,options)
							for Peak=1:NumPeaks
									if TrialParameters(2*Peak)<MINWIDTH,
											TrialParameters(2*Peak)=MINWIDTH
									
							
					case 15
							TrialParameters=fminsearch(@(lambda)(fitBWF(lambda,xx,yy,extra)),newstart,options)
							for Peak=1:NumPeaks
									if TrialParameters(2*Peak)<MINWIDTH,
											TrialParameters(2*Peak)=MINWIDTH
									
							
					case 16
							fixedstart=[]
							for pc=1:NumPeaks,
									fixedstart(pc)=(max(xx)-min(xx))./(NumPeaks+1)
							
							TrialParameters=fminsearch(@(lambda)(FitFPGaussian(lambda,xx,yy)),fixedstart,options)
							for Peak=1:NumPeaks
									if TrialParameters(Peak)<MINWIDTH,
											TrialParameters(Peak)=MINWIDTH
									
							
					case 17
							fixedstart=[]
							for pc=1:NumPeaks,
									fixedstart(pc)=(max(xx)-min(xx))./(NumPeaks+1)
							
							TrialParameters=fminsearch(@(lambda)(FitFPLorentzian(lambda,xx,yy)),fixedstart,options)
							for Peak=1:NumPeaks
									if TrialParameters(Peak)<MINWIDTH,
											TrialParameters(Peak)=MINWIDTH
									
							
					case 18
							zxx=[zeros(size(xx)) xx zeros(size(xx)) ]
							zyy=[zeros(size(yy)) yy zeros(size(yy)) ]
							TrialParameters=fminsearch(@(lambda)(fitexplorentzian(lambda,zxx,zyy,-extra)),newstart,options)
					case 19
							TrialParameters=fminsearch(@(lambda)(alphadef(lambda,xx,yy)),newstart,options)
							for Peak=1:NumPeaks
									if TrialParameters(2*Peak)<MINWIDTH,
											TrialParameters(2*Peak)=MINWIDTH
									
							
					case 20
							TrialParameters=fminsearch(@(lambda)(fitvoigt(lambda,xx,yy,extra)),newstart,options)
							for Peak=1:NumPeaks
									if TrialParameters(2*Peak)<MINWIDTH,
											TrialParameters(2*Peak)=MINWIDTH
									
							
					case 21
							TrialParameters=fminsearch(@(lambda)(fittriangular(lambda,xx,yy)),newstart,options)
							for Peak=1:NumPeaks
									if TrialParameters(2*Peak)<MINWIDTH,
											TrialParameters(2*Peak)=MINWIDTH
									
							
					case 22
							TrialParameters=fminsearch(@(lambda)(fitmultiple(lambda,xx,yy,shapesvector,extra)),newstart,options)
							for Peak=1:NumPeaks
									if TrialParameters(2*Peak)<MINWIDTH(Peak),
											TrialParameters(2*Peak)=MINWIDTH(Peak)
									
							
					case 23
							TrialParameters=fminsearch(@(lambda)(fitdownsigmoid(lambda,xx,yy)),newstar,optionst)
							for Peak=1:NumPeaks
									if TrialParameters(2*Peak)<MINWIDTH,
											TrialParameters(2*Peak)=MINWIDTH
									
							
					case 24
							TrialParameters=fminsearch(@(lambda)(fitnbinpdf(lambda,xx,yy)),newstart,options)
							for Peak=1:NumPeaks
									if TrialParameters(2*Peak)<MINWIDTH,
											TrialParameters(2*Peak)=MINWIDTH
									
							
					case 25
							TrialParameters=fminsearch(@(lambda)(fitlognpdf(lambda,xx,yy)),newstart,options)
							for Peak=1:NumPeaks
									if TrialParameters(2*Peak)<MINWIDTH,
											TrialParameters(2*Peak)=MINWIDTH
									
							
					case 26
							TrialParameters=fminsearch(@(lambda)(fitlinslope(lambda,xx,yy)),polyfit(xx,yy,1),options)
					coeff=TrialParameters
					case 27
							TrialParameters=fminsearch(@(lambda)(fitd1gauss(lambda,xx,yy)),newstart,options)
							for Peak=1:NumPeaks
									if TrialParameters(2*Peak)<MINWIDTH,
											TrialParameters(2*Peak)=MINWIDTH
									
							
					case 28
							TrialParameters=fitpolynomial(xx,yy,extra)
					case 29
							TrialParameters=fminsearch(@(lambda)(fitsegmented(lambda,xx,yy)),newstart,options)
					case 30
							TrialParameters=fminsearch(@(lambda)(fitvoigtv(lambda,xx,yy)),newstart)
					case 31
							zxx=[zeros(size(xx)) xx zeros(size(xx)) ]
							zyy=[zeros(size(yy)) yy zeros(size(yy)) ]
							TrialParameters=fminsearch(@(lambda)(fitexpgaussianv(lambda,zxx,zyy)),newstart)
					case 32
							TrialParameters=fminsearch(@(lambda)(fitpearsonv(lambda,xx,yy)),newstart)
					case 33
							TrialParameters=fminsearch(@(lambda)(fitGLv(lambda,xx,yy)),newstart)
					otherwise
			 # switch peakshape
			
		for peaks=1:NumPeaks,
			 peakindex=2*peaks-1
			 newstart(peakindex)=start(peakindex)-xoffset


			# Construct model from Trial parameters
			A=zeros(NumPeaks,n)
			for m=1:NumPeaks,
					switch peakshape(1)
							case 1
									A(m,:)=gaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m))
							case 2
									A(m,:)=lorentzian(xx,TrialParameters(2*m-1),TrialParameters(2*m))
							case 3
									A(m,:)=logistic(xx,TrialParameters(2*m-1),TrialParameters(2*m))
							case 4
									A(m,:)=pearson(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra)
							case 5
									A(m,:)=expgaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m),-extra)'
							case 6
									A(m,:)=gaussian(xx,TrialParameters(m),TrialParameters(NumPeaks+1))
							case 7
									A(m,:)=lorentzian(xx,TrialParameters(m),TrialParameters(NumPeaks+1))
							case 8
									A(m,:)=expgaussian(xx,TrialParameters(m),TrialParameters(NumPeaks+1),-extra)'
							case 9
									A(m,:)=exppulse(xx,TrialParameters(2*m-1),TrialParameters(2*m))
							case 10
									A(m,:)=upsigmoid(xx,TrialParameters(2*m-1),TrialParameters(2*m))
							case 11
									A(m,:)=gaussian(xx,TrialParameters(m),FIXEDPARAMETERS)
							case 12
									A(m,:)=lorentzian(xx,TrialParameters(m),FIXEDPARAMETERS)
							case 13
									A(m,:)=GL(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra)
							case 14
									A(m,:)=BiGaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra)
							case 15
									A(m,:)=BWF(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra)
							case 16
									A(m,:)=gaussian(xx,FIXEDPOSITIONS(m),TrialParameters(m))
							case 17
									A(m,:)=lorentzian(xx,FIXEDPOSITIONS(m),TrialParameters(m))
							case 18
									A(m,:)=explorentzian(xx,TrialParameters(2*m-1),TrialParameters(2*m),-extra)'
							case 19
									A(m,:)=alphadef(xx,TrialParameters(2*m-1),TrialParameters(2*m))
							case 20
									A(m,:)=voigt(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra)
							case 21
									A(m,:)=triangular(xx,TrialParameters(2*m-1),TrialParameters(2*m))
							case 22
									A(m,:)=peakdef(shapesvector(m),xx,TrialParameters(2*m-1),TrialParameters(2*m),extra(m))
							case 23
									A(m,:)=downsigmoid(xx,TrialParameters(2*m-1),TrialParameters(2*m))      
							case 24
									A(m,:)=nbinpdf(xx,TrialParameters(2*m-1),TrialParameters(2*m))
							case 25
									A(m,:)=lognormal(xx,TrialParameters(2*m-1),TrialParameters(2*m))
							case 26
									A(m,:)=linslope(xx,TrialParameters(2*m-1),TrialParameters(2*m))
							case 27
									A(m,:)=d1gauss(xx,TrialParameters(2*m-1),TrialParameters(2*m))       
							case 28
									A(m,:)=polynomial(xx,TrialParameters(2*m-1),TrialParameters(2*m))       
							case 29
									A(m,:)=segmented(xx,yy,PEAKHEIGHTS)
							case 30
									A(m,:)=voigt(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m))        
							case 31
									A(m,:)=expgaussian(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m))        
							case 32
									A(m,:)=pearson(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m))        
							case 33
									A(m,:)=GL(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m))
					 # switch
			 # for
			
			# Multiplies each row by the corresponding amplitude and adds them up
			if peakshape(1)==29, # Segmented linear
					model=segmented(xx,yy,PEAKHEIGHTS)
					TrialParameters=coeff
					Heights=ones(size(coeff))
			else:
					if AUTOZERO==3,
							baseline=PEAKHEIGHTS(1)
							Heights=PEAKHEIGHTS(2:1+NumPeaks)
							model=Heights'*A+baseline
					else:
							model=PEAKHEIGHTS'*A
							Heights=PEAKHEIGHTS
							baseline=0
					
			
			
			# Compare trial model to data segment and compute the fit error
			MeanFitError=100*norm(yy-model)./(sqrt(n)*max(yy))
			# Take only the single fit that has the lowest MeanFitError
			if MeanFitError<LowestError,
					if min(Heights)>=-BIPOLAR*10^100,  # Consider only fits with positive peak heights
							LowestError=MeanFitError  # Assign LowestError to the lowest MeanFitError
							FitParameters=TrialParameters  # Assign FitParameters to the fit with the lowest MeanFitError
							height=Heights # Assign height to the PEAKHEIGHTS with the lowest MeanFitError
					 # if min(PEAKHEIGHTS)>0
			 # if MeanFitError<LowestError
		# for k (NumTrials)
			Rsquared=1-(norm(yy-bestmodel)./norm(yy-mean(yy)))
			SStot=sum((yy-mean(yy)).^2)
			SSres=sum((yy-bestmodel).^2)
			Rsquared=1-(SSres./SStot)
			GOF=[LowestError Rsquared]
		for m=1:NumPeaks,
			area(m)=trapz(xx+xoffset,height(m)*A(m,:)) # Compute the area of each component peak using trapezoidal method


		# Put results into a matrix FitResults, one row for each peak, showing peak index number,
		# position, amplitude, and width.
		FitResults=zeros(NumPeaks,6)
		switch peakshape(1),
			case {6,7,8}, # equal-width peak models only
					for m=1:NumPeaks,
							if m==1,
									FitResults=[[round(m) FitParameters(m)+xoffset height(m) abs(FitParameters(NumPeaks+1)) area(m)]]
							else:
									FitResults=[FitResults  [round(m) FitParameters(m)+xoffset height(m) abs(FitParameters(NumPeaks+1)) area(m)]]
							
					
			case {11,12}, # Fixed-width shapes only
					for m=1:NumPeaks,
							if m==1,
									FitResults=[[round(m) FitParameters(m)+xoffset height(m) FIXEDPARAMETERS area(m)]]
							else:
									FitResults=[FitResults  [round(m) FitParameters(m)+xoffset height(m) FIXEDPARAMETERS area(m)]]
							
					
			case {16,17}, # Fixed-position shapes only
					for m=1:NumPeaks,
							if m==1,
									FitResults=[round(m) FIXEDPOSITIONS(m) height(m) FitParameters(m) area(m)]
							else:
									FitResults=[FitResults  [round(m) FIXEDPOSITIONS(m) height(m) FitParameters(m) area(m)]]
							
					
			case 28,   # Simple polynomial fit
					FitResults=PEAKHEIGHTS
			case 29, # Segmented linear fit
					FitResults=PEAKHEIGHTS
			case {30,31,32,33} # Special case of shapes with 3 iterated variables
					for m=1:NumPeaks,
							if m==1,
									FitResults=[round(m) FitParameters(3*m-2) height(m) abs(FitParameters(3*m-1)) area(m) FitParameters(3*m)]
							else:
									FitResults=[FitResults  [round(m) FitParameters(3*m-2) height(m) abs(FitParameters(3*m-1)) area(m) FitParameters(3*m)]]
							
					
			otherwise # Normal shapes with 2 iterated variables
					for m=1:NumPeaks,
							if m==1,
									FitResults=[round(m) FitParameters(2*m-1)+xoffset height(m) abs(FitParameters(2*m)) area(m)]
							else:
									FitResults=[FitResults  [round(m) FitParameters(2*m-1)+xoffset height(m) abs(FitParameters(2*m)) area(m)]]
							 # if m==1
					 # for m=1:NumPeaks,
		# switch peakshape(1)
		# ----------------------------------------------------------------------
		def start=calcstart(xx,NumPeaks,xoffset)
		n=max(xx)-min(xx)
		start=[]
		startpos=[n/(NumPeaks+1):n/(NumPeaks+1):n-(n/(NumPeaks+1))]+min(xx)
		for marker=1:NumPeaks,
				markx=startpos(marker)+ xoffset
				start=[start markx n/ (3.*NumPeaks)]
		 # for marker
		# ----------------------------------------------------------------------
		def err = fitgaussian(lambda,t,y)
		# Fitting def for a Gaussian band signal.
		global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
		numpeaks=round(length(lambda)/2)
		A = zeros(length(t),numpeaks)
		for j = 1:numpeaks,
		#    if lambda(2*j)<MINWIDTH,lambda(2*j)=MINWIDTH
			A(:,j) = gaussian(t,lambda(2*j-1),lambda(2*j))'

		if AUTOZERO==3,A=[ones(size(y))' A]
		if BIPOLAR,PEAKHEIGHTS=A\y'else: PEAKHEIGHTS=abs(A\y')
		z = A*PEAKHEIGHTS
		if LOGPLOT,
			err = norm(log10(z)-log10(y)')
		else:
			err = norm(z-y')

		# ----------------------------------------------------------------------
		def err = fitewgaussian(lambda,t,y)
		# Fitting def for a Gaussian band signal with equal peak widths.
		global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
		numpeaks=round(length(lambda)-1)
		A = zeros(length(t),numpeaks)
		for j = 1:numpeaks,
			A(:,j) = gaussian(t,lambda(j),lambda(numpeaks+1))'

		if AUTOZERO==3,A=[ones(size(y))' A]
		if BIPOLAR,PEAKHEIGHTS=A\y'else: PEAKHEIGHTS=abs(A\y')
		z = A*PEAKHEIGHTS
		if LOGPLOT,
			err = norm(log10(z)-log10(y)')
		else:
			err = norm(z-y')

		# ----------------------------------------------------------------------
		def err = FitFWGaussian(lambda,t,y)
		#	Fitting def for a fixed width Gaussian
		global PEAKHEIGHTS AUTOZERO FIXEDPARAMETERS BIPOLAR LOGPLOT
		numpeaks=round(length(lambda))
		A = zeros(length(t),numpeaks)
		for j = 1:numpeaks,
			A(:,j) = gaussian(t,lambda(j),FIXEDPARAMETERS)'

		if AUTOZERO==3,A=[ones(size(y))' A]
		if BIPOLAR,PEAKHEIGHTS=A\y'else: PEAKHEIGHTS=abs(A\y')
		z = A*PEAKHEIGHTS
		if LOGPLOT,
			err = norm(log10(z)-log10(y)')
		else:
			err = norm(z-y')

		# ----------------------------------------------------------------------
		def err = FitFPGaussian(lambda,t,y)
		#	Fitting def for fixed-position Gaussians
		global PEAKHEIGHTS AUTOZERO FIXEDPARAMETERS BIPOLAR LOGPLOT
		numpeaks=round(length(lambda))
		A = zeros(length(t),numpeaks)
		for j = 1:numpeaks,
			A(:,j) = gaussian(t,FIXEDPARAMETERS(j), lambda(j))'

		if AUTOZERO==3,A=[ones(size(y))' A]
		if BIPOLAR,PEAKHEIGHTS=A\y'else: PEAKHEIGHTS=abs(A\y')
		z = A*PEAKHEIGHTS
		if LOGPLOT,
			err = norm(log10(z)-log10(y)')
		else:
			err = norm(z-y')

		# ----------------------------------------------------------------------
		def err = FitFPLorentzian(lambda,t,y)
		#	Fitting def for fixed-position Lorentzians
		global PEAKHEIGHTS AUTOZERO FIXEDPARAMETERS BIPOLAR
		numpeaks=round(length(lambda))
		A = zeros(length(t),numpeaks)
		for j = 1:numpeaks,
			A(:,j) = lorentzian(t,FIXEDPARAMETERS(j), lambda(j))'

		if AUTOZERO==3,A=[ones(size(y))' A]
		if BIPOLAR,PEAKHEIGHTS=A\y'else: PEAKHEIGHTS=abs(A\y')
		z = A*PEAKHEIGHTS
		err = norm(z-y')
		# ----------------------------------------------------------------------
		def err = FitFWLorentzian(lambda,t,y)
		#	Fitting def for fixed width Lorentzian
		global PEAKHEIGHTS AUTOZERO FIXEDPARAMETERS BIPOLAR LOGPLOT
		numpeaks=round(length(lambda))
		A = zeros(length(t),numpeaks)
		for j = 1:numpeaks,
			A(:,j) = lorentzian(t,lambda(j),FIXEDPARAMETERS)'

		if AUTOZERO==3,A=[ones(size(y))' A]
		if BIPOLAR,PEAKHEIGHTS=A\y'else: PEAKHEIGHTS=abs(A\y')
		z = A*PEAKHEIGHTS
		if LOGPLOT,
			err = norm(log10(z)-log10(y)')
		else:
			err = norm(z-y')

		# ----------------------------------------------------------------------
		def err = fitewlorentzian(lambda,t,y)
		# Fitting def for a Lorentzian band signal with equal peak widths.
		global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
		numpeaks=round(length(lambda)-1)
		A = zeros(length(t),numpeaks)
		for j = 1:numpeaks,
			A(:,j) = lorentzian(t,lambda(j),lambda(numpeaks+1))'

		if AUTOZERO==3,A=[ones(size(y))' A]
		if BIPOLAR,PEAKHEIGHTS=A\y'else: PEAKHEIGHTS=abs(A\y')
		z = A*PEAKHEIGHTS
		if LOGPLOT,
			err = norm(log10(z)-log10(y)')
		else:
			err = norm(z-y')

		# ----------------------------------------------------------------------
		def g = gaussian(x,pos,wid)
		#  gaussian(X,pos,wid) = gaussian peak centered on pos, half-width=wid
		#  X may be scalar, vector, or matrix, pos and wid both scalar
		# Examples: gaussian([0 1 2],1,2) gives result [0.5000    1.0000    0.5000]
		# plot(gaussian([1:100],50,20)) displays gaussian band centered at 50 with width 20.
		g = exp(-((x-pos)./(0.6005615.*wid)).^2)
		# ----------------------------------------------------------------------
		def err = fitlorentzian(lambda,t,y)
		#	Fitting def for single lorentzian, lambda(1)=position, lambda(2)=width
		#	Fitgauss assumes a lorentzian def 
		global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
		A = zeros(length(t),round(length(lambda)/2))
		for j = 1:length(lambda)/2,
			A(:,j) = lorentzian(t,lambda(2*j-1),lambda(2*j))'

		if AUTOZERO==3,A=[ones(size(y))' A]
		if BIPOLAR,PEAKHEIGHTS=A\y'else: PEAKHEIGHTS=abs(A\y')
		z = A*PEAKHEIGHTS
		if LOGPLOT,
			err = norm(log10(z)-log10(y)')
		else:
			err = norm(z-y')

		# ----------------------------------------------------------------------
		def g = lorentzian(x,position,width)
		# lorentzian(x,position,width) Lorentzian def.
		# where x may be scalar, vector, or matrix
		# position and width scalar
		# T. C. O'Haver, 1988
		# Example: lorentzian([1 2 3],2,2) gives result [0.5 1 0.5]
		g=ones(size(x))./(1+((x-position)./(0.5.*width)).^2)
		# ----------------------------------------------------------------------
		def err = fitlogistic(lambda,t,y)
		#	Fitting def for logistic, lambda(1)=position, lambda(2)=width
		#	between the data and the values computed by the current
		#	def of lambda.  Fitlogistic assumes a logistic def 
		#  T. C. O'Haver, May 2006
		global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
		A = zeros(length(t),round(length(lambda)/2))
		for j = 1:length(lambda)/2,
			A(:,j) = logistic(t,lambda(2*j-1),lambda(2*j))'

		if AUTOZERO==3,A=[ones(size(y))' A]
		if BIPOLAR,PEAKHEIGHTS=A\y'else: PEAKHEIGHTS=abs(A\y')
		z = A*PEAKHEIGHTS
		if LOGPLOT,
			err = norm(log10(z)-log10(y)')
		else:
			err = norm(z-y')

		# ----------------------------------------------------------------------
		def g = logistic(x,pos,wid)
		# logistic def.  pos=position wid=half-width (both scalar)
		# logistic(x,pos,wid), where x may be scalar, vector, or matrix
		# pos=position wid=half-width (both scalar)
		# T. C. O'Haver, 1991 
		n = exp(-((x-pos)/(.477.*wid)) .^2)
		g = (2.*n)./(1+n)
		# ----------------------------------------------------------------------
		def err = fittriangular(lambda,t,y)
		#	Fitting def for triangular, lambda(1)=position, lambda(2)=width
		#	between the data and the values computed by the current
		#	def of lambda.  Fittriangular assumes a triangular def 
		#  T. C. O'Haver, May 2006
		global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
		A = zeros(length(t),round(length(lambda)/2))
		for j = 1:length(lambda)/2,
			A(:,j) = triangular(t,lambda(2*j-1),lambda(2*j))'

		if AUTOZERO==3,A=[ones(size(y))' A]
		if BIPOLAR,PEAKHEIGHTS=A\y'else: PEAKHEIGHTS=abs(A\y')
		z = A*PEAKHEIGHTS
		if LOGPLOT,
			err = norm(log10(z)-log10(y)')
		else:
			err = norm(z-y')

		# ----------------------------------------------------------------------
		def g = triangular(x,pos,wid)
		#triangle def.  pos=position wid=half-width (both scalar)
		#trianglar(x,pos,wid), where x may be scalar or vector,
		#pos=position wid=half-width (both scalar)
		# T. C. O'Haver, 1991
		# Example
		# x=[0:.1:10]plot(x,trianglar(x,5.5,2.3),'.')
		g=1-(1./wid) .*abs(x-pos)
		for i=1:length(x),  
		if g(i)<0,g(i)=0

		# ----------------------------------------------------------------------
		def err = fitpearson(lambda,t,y,shapeconstant)
		#   Fitting defs for a Pearson 7 band signal.
		# T. C. O'Haver (toh@umd.edu),   Version 1.3, October 23, 2006.
		global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
		A = zeros(length(t),round(length(lambda)/2))
		for j = 1:length(lambda)/2,
			A(:,j) = pearson(t,lambda(2*j-1),lambda(2*j),shapeconstant)'

		if AUTOZERO==3,A=[ones(size(y))' A]
		if BIPOLAR,PEAKHEIGHTS=A\y'else: PEAKHEIGHTS=abs(A\y')
		z = A*PEAKHEIGHTS
		if LOGPLOT,
			err = norm(log10(z)-log10(y)')
		else:
			err = norm(z-y')

		# ----------------------------------------------------------------------
		def err = fitpearsonv(lambda,t,y)
		# Fitting defs for pearson def with indepently variable
		# percent Gaussian
		# T. C. O'Haver (toh@umd.edu), 2015.
		global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
		A = zeros(length(t),round(length(lambda)/3))
		for j = 1:length(lambda)/3,
			A(:,j) = pearson(t,lambda(3*j-2),lambda(3*j-1),lambda(3*j))'

		if AUTOZERO==3,A=[ones(size(y))' A]
		if BIPOLAR,PEAKHEIGHTS=A\y'else: PEAKHEIGHTS=abs(A\y')
		z = A*PEAKHEIGHTS
		if LOGPLOT,
			err = norm(log10(z)-log10(y)')
		else:
			err = norm(z-y')

		# ----------------------------------------------------------------------
		def g = pearson(x,pos,wid,m)
		# Pearson VII def. 
		# g = pearson7(x,pos,wid,m) where x may be scalar, vector, or matrix
		# pos=position wid=half-width (both scalar)
		# m=some number
		#  T. C. O'Haver, 1990  
		g=ones(size(x))./(1+((x-pos)./((0.5.^(2/m)).*wid)).^2).^m
		# ----------------------------------------------------------------------
		def err = fitexpgaussian(lambda,t,y,timeconstant)
		#   Fitting defs for a exponentially-broadened Gaussian band signal.
		#  T. C. O'Haver, October 23, 2006.
		global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
		A = zeros(length(t),round(length(lambda)/2))
		for j = 1:length(lambda)/2,
			A(:,j) = expgaussian(t,lambda(2*j-1),lambda(2*j),timeconstant)

		if AUTOZERO==3,A=[ones(size(y))' A]
		if BIPOLAR,PEAKHEIGHTS=A\y'else: PEAKHEIGHTS=abs(A\y')
		z = A*PEAKHEIGHTS
		if LOGPLOT,
			err = norm(log10(z)-log10(y)')
		else:
			err = norm(z-y')

		# ----------------------------------------------------------------------
		def err = fitexplorentzian(lambda,t,y,timeconstant)
		#   Fitting defs for a exponentially-broadened lorentzian band signal.
		#  T. C. O'Haver, 2013.
		global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
		A = zeros(length(t),round(length(lambda)/2))
		for j = 1:length(lambda)/2,
			A(:,j) = explorentzian(t,lambda(2*j-1),lambda(2*j),timeconstant)

		if AUTOZERO==3,A=[ones(size(y))' A]
		if BIPOLAR,PEAKHEIGHTS=A\y'else: PEAKHEIGHTS=abs(A\y')
		z = A*PEAKHEIGHTS
		if LOGPLOT,
			err = norm(log10(z)-log10(y)')
		else:
			err = norm(z-y')

		# ----------------------------------------------------------------------
		def err = fitexpewgaussian(lambda,t,y,timeconstant)
		# Fitting def for exponentially-broadened Gaussian bands with equal peak widths.
		global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
		numpeaks=round(length(lambda)-1)
		A = zeros(length(t),numpeaks)
		for j = 1:numpeaks,
			A(:,j) = expgaussian(t,lambda(j),lambda(numpeaks+1),timeconstant)

		if AUTOZERO==3,A=[ones(size(y))' A]
		if BIPOLAR,PEAKHEIGHTS=A\y'else: PEAKHEIGHTS=abs(A\y')
		z = A*PEAKHEIGHTS
		if LOGPLOT,
			err = norm(log10(z)-log10(y)')
		else:
			err = norm(z-y')

		# ----------------------------------------------------------------------
		def err = fitexpgaussianv(lambda,t,y)
		# Fitting defs for  exponentially-broadened Gaussians with
		# indepently variable time constants
		# T. C. O'Haver (toh@umd.edu), 2015.
		global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
		A = zeros(length(t),round(length(lambda)/3))
		for j = 1:length(lambda)/3,
			A(:,j) = expgaussian(t,lambda(3*j-2),lambda(3*j-1),-lambda(3*j))'

		if AUTOZERO==3,A=[ones(size(y))' A]
		if BIPOLAR,PEAKHEIGHTS=A\y'else: PEAKHEIGHTS=abs(A\y')
		z = A*PEAKHEIGHTS
		if LOGPLOT,
			err = norm(log10(z)-log10(y)')
		else:
			err = norm(z-y')

		# ----------------------------------------------------------------------
		def g = expgaussian(x,pos,wid,timeconstant)
		#  Exponentially-broadened gaussian(x,pos,wid) = gaussian peak centered on pos, half-width=wid
		#  x may be scalar, vector, or matrix, pos and wid both scalar
		#  T. C. O'Haver, 2006
		g = exp(-((x-pos)./(0.6005615.*wid)) .^2)
		g = ExpBroaden(g',timeconstant)
		# ----------------------------------------------------------------------
		def g = explorentzian(x,pos,wid,timeconstant)
		#  Exponentially-broadened lorentzian(x,pos,wid) = lorentzian peak centered on pos, half-width=wid
		#  x may be scalar, vector, or matrix, pos and wid both scalar
		#  T. C. O'Haver, 2013
		g = ones(size(x))./(1+((x-pos)./(0.5.*wid)).^2)
		g = ExpBroaden(g',timeconstant)
		# ----------------------------------------------------------------------
		def yb = ExpBroaden(y,t)
		# ExpBroaden(y,t) zero pads y and convolutes result by an exponential decay
		# of time constant t by multiplying Fourier transforms and inverse
		# transforming the result.
		hly=round(length(y)./2)
		ey=[y(1).*ones(1,hly)'yy(length(y)).*ones(1,hly)']
		# figure(2)plot(ey)figure(1)
		fy=fft(ey)
		a=exp(-(1:length(fy))./t)
		fa=fft(a)
		fy1=fy.*fa'
		ybz=real(ifft(fy1))./sum(a)
		yb=ybz(hly+2:length(ybz)-hly+1)
		# ----------------------------------------------------------------------
		def err = fitexppulse(tau,x,y)
		# Iterative fit of the sum of exponential pulses
		# of the form Height.*exp(-tau1.*x).*(1-exp(-tau2.*x)))
		global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
		A = zeros(length(x),round(length(tau)/2))
		for j = 1:length(tau)/2,
			A(:,j) = exppulse(x,tau(2*j-1),tau(2*j))

		if AUTOZERO==3,A=[ones(size(y))' A]
		if BIPOLAR,PEAKHEIGHTS=A\y'else: PEAKHEIGHTS=abs(A\y')
		z = A*PEAKHEIGHTS
		if LOGPLOT,
			err = norm(log10(z)-log10(y)')
		else:
			err = norm(z-y')

		# ----------------------------------------------------------------------
		def g = exppulse(x,t1,t2)
		# Exponential pulse of the form 
		# g = (x-spoint)./pos.*exp(1-(x-spoint)./pos)
		e=(x-t1)./t2
		p = 4*exp(-e).*(1-exp(-e))
		p=p .* (p>0)
		g = p'
		# ----------------------------------------------------------------------
		def err = fitalphadef(tau,x,y)
		# Iterative fit of the sum of alpha funciton
		# of the form Height.*exp(-tau1.*x).*(1-exp(-tau2.*x)))
		global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
		A = zeros(length(x),round(length(tau)/2))
		for j = 1:length(tau)/2,
			A(:,j) = alphadef(x,tau(2*j-1),tau(2*j))

		if AUTOZERO==3,A=[ones(size(y))' A]
		if BIPOLAR,PEAKHEIGHTS=A\y'else: PEAKHEIGHTS=abs(A\y')
		z = A*PEAKHEIGHTS
		if LOGPLOT,
			err = norm(log10(z)-log10(y)')
		else:
			err = norm(z-y')

		# ----------------------------------------------------------------------
		def g = alphadef(x,pos,spoint)
		# alpha def.  pos=position wid=half-width (both scalar)
		# alphadef(x,pos,wid), where x may be scalar, vector, or matrix
		# pos=position wid=half-width (both scalar)
		# Taekyung Kwon, July 2013  
		g = (x-spoint)./pos.*exp(1-(x-spoint)./pos)
		for m=1:length(x)if g(m)<0g(m)=0
		# ----------------------------------------------------------------------
		def err = fitdownsigmoid(tau,x,y)
		# Fitting def for iterative fit to the sum of
		# downward moving sigmiods 
		global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
		A = zeros(length(x),round(length(tau)/2))
		for j = 1:length(tau)/2,
			A(:,j) = downsigmoid(x,tau(2*j-1),tau(2*j))

		if AUTOZERO==3,A=[ones(size(y))' A]
		if BIPOLAR,PEAKHEIGHTS=A\y'else: PEAKHEIGHTS=abs(A\y')
		z = A*PEAKHEIGHTS
		if LOGPLOT,
			err = norm(log10(z)-log10(y)')
		else:
			err = norm(z-y')

		# ----------------------------------------------------------------------
		def err = fitupsigmoid(tau,x,y)
		# Fitting def for iterative fit to the sum of
		# upwards moving sigmiods
		global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
		A = zeros(length(x),round(length(tau)/2))
		for j = 1:length(tau)/2,
			A(:,j) = upsigmoid(x,tau(2*j-1),tau(2*j))

		if AUTOZERO==3,A=[ones(size(y))' A]
		if BIPOLAR,PEAKHEIGHTS=A\y'else: PEAKHEIGHTS=abs(A\y')
		z = A*PEAKHEIGHTS
		if LOGPLOT,
			err = norm(log10(z)-log10(y)')
		else:
			err = norm(z-y')

		# ----------------------------------------------------------------------
		def g=downsigmoid(x,t1,t2)
		# down step sigmoid
		g=.5-.5*erf(real((x-t1)/sqrt(2*t2)))
		# ----------------------------------------------------------------------
		def g=upsigmoid(x,t1,t2)
		# up step sigmoid
		g=1/2 + 1/2* erf(real((x-t1)/sqrt(2*t2))) 
		# ----------------------------------------------------------------------
		def err = fitGL(lambda,t,y,shapeconstant)
		#   Fitting defs for Gaussian/Lorentzian bl.
		# T. C. O'Haver (toh@umd.edu), 2012.
		global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
		A = zeros(length(t),round(length(lambda)/2))
		for j = 1:length(lambda)/2,
			A(:,j) = GL(t,lambda(2*j-1),lambda(2*j),shapeconstant)'

		if AUTOZERO==3,A=[ones(size(y))' A]
		if BIPOLAR,PEAKHEIGHTS=A\y'else: PEAKHEIGHTS=abs(A\y')
		z = A*PEAKHEIGHTS
		if LOGPLOT,
			err = norm(log10(z)-log10(y)')
		else:
			err = norm(z-y')

		# ----------------------------------------------------------------------
		def err = fitGLv(lambda,t,y)
		# Fitting defs for Gaussian/Lorentzian bl def with
		# indepently variable percent Gaussian
		# T. C. O'Haver (toh@umd.edu), 2015.
		global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
		A = zeros(length(t),round(length(lambda)/3))
		for j = 1:length(lambda)/3,
			A(:,j) = GL(t,lambda(3*j-2),lambda(3*j-1),lambda(3*j))'

		if AUTOZERO==3,A=[ones(size(y))' A]
		if BIPOLAR,PEAKHEIGHTS=A\y'else: PEAKHEIGHTS=abs(A\y')
		z = A*PEAKHEIGHTS
		if LOGPLOT,
			err = norm(log10(z)-log10(y)')
		else:
			err = norm(z-y')

		# ----------------------------------------------------------------------
		def g = GL(x,pos,wid,m)
		# Gaussian/Lorentzian bl. m = percent Gaussian character
		# pos=position wid=half-width
		# m = percent Gaussian character.
		#  T. C. O'Haver, 2012
		# sizex=size(x)
		# sizepos=size(pos)
		# sizewid=size(wid)
		# sizem=size(m)
		g=2.*((m/100).*gaussian(x,pos,wid)+(1-(m(1)/100)).*lorentzian(x,pos,wid))/2
		# ----------------------------------------------------------------------
		def err = fitvoigt(lambda,t,y,shapeconstant)
		# Fitting defs for Voigt profile def
		# T. C. O'Haver (toh@umd.edu), 2013.
		global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
		A = zeros(length(t),round(length(lambda)/2))
		for j = 1:length(lambda)/2,
			A(:,j) = voigt(t,lambda(2*j-1),lambda(2*j),shapeconstant)'

		if AUTOZERO==3,A=[ones(size(y))' A]
		if BIPOLAR,PEAKHEIGHTS=A\y'else: PEAKHEIGHTS=abs(A\y')
		z = A*PEAKHEIGHTS
		if LOGPLOT,
			err = norm(log10(z)-log10(y)')
		else:
			err = norm(z-y')

		# ----------------------------------------------------------------------
		def err = fitvoigtv(lambda,t,y)
		# Fitting defs for Voigt profile def with indepently variable
		# alphas
		# T. C. O'Haver (toh@umd.edu), 2015.
		global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
		A = zeros(length(t),round(length(lambda)/3))
		for j = 1:length(lambda)/3,
			A(:,j) = voigt(t,lambda(3*j-2),lambda(3*j-1),lambda(3*j))'

		if AUTOZERO==3,A=[ones(size(y))' A]
		if BIPOLAR,PEAKHEIGHTS=A\y'else: PEAKHEIGHTS=abs(A\y')
		z = A*PEAKHEIGHTS
		if LOGPLOT,
			err = norm(log10(z)-log10(y)')
		else:
			err = norm(z-y')

		# ----------------------------------------------------------------------
		def g=voigt(xx,pos,gD,alpha)
		# Voigt profile def. xx is the indepent variable (energy,
		# wavelength, etc), gD is the Doppler (Gaussian) width, and alpha is the
		# shape constant (ratio of the Lorentzian width gL to the Doppler width gD.
		# Based on Chong Tao's "Voigt lineshape spectrum simulation", 
		# File ID: #26707
		% alpha=alpha
		gL=alpha.*gD
		gV = 0.5346*gL + sqrt(0.2166*gL.^2 + gD.^2)
		x = gL/gV
		y = abs(xx-pos)/gV
		g = 1/(2*gV*(1.065 + 0.447*x + 0.058*x^2))*((1-x)*exp(-0.693.*y.^2) + (x./(1+y.^2)) + 0.016*(1-x)*x*(exp(-0.0841.*y.^2.25)-1./(1 + 0.021.*y.^2.25)))
		g=g./max(g)
		% ----------------------------------------------------------------------
		def err = fitBiGaussian(lambda,t,y,shapeconstant)
		%   Fitting defs for BiGaussian.
		% T. C. O'Haver (toh@umd.edu),  2012.
		global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
		A = zeros(length(t),round(length(lambda)/2))
		for j = 1:length(lambda)/2,
			A(:,j) = BiGaussian(t,lambda(2*j-1),lambda(2*j),shapeconstant)'

		if AUTOZERO==3,A=[ones(size(y))' A] 
		if BIPOLAR,PEAKHEIGHTS=A\y'else: PEAKHEIGHTS=abs(A\y')
		z = A*PEAKHEIGHTS
		if LOGPLOT,
			err = norm(log10(z)-log10(y)')
		else:
			err = norm(z-y')

		% ----------------------------------------------------------------------
		def g = BiGaussian(x,pos,wid,m)
		% BiGaussian (different widths on leading edge and trailing edge).
		% pos=position wid=width 
		% m determines shape symmetrical if m=50.
		%  T. C. O'Haver, 2012
		lx=length(x)
		hx=val2ind(x,pos)
		g(1:hx)=gaussian(x(1:hx),pos,wid*(m/100))
		g(hx+1:lx)=gaussian(x(hx+1:lx),pos,wid*(1-m/100))
		% ----------------------------------------------------------------------
		def err = fitBWF(lambda,t,y,shapeconstant)
		%   Fitting def for Breit-Wigner-Fano.
		% T. C. O'Haver (toh@umd.edu),  2014.
		global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
		A = zeros(length(t),round(length(lambda)/2))
		for j = 1:length(lambda)/2,
			A(:,j) = BWF(t,lambda(2*j-1),lambda(2*j),shapeconstant)'

		if AUTOZERO==3,A=[ones(size(y))' A]
		if BIPOLAR,PEAKHEIGHTS=A\y'else: PEAKHEIGHTS=abs(A\y')
		z = A*PEAKHEIGHTS
		if LOGPLOT,
			err = norm(log10(z)-log10(y)')
		else:
			err = norm(z-y')

		% ----------------------------------------------------------------------
		def g = BWF(x,pos,wid,m)
		% BWF (Breit-Wigner-Fano) http://en.wikipedia.org/wiki/Fano_resonance
		% pos=position wid=width m=Fano factor
		%  T. C. O'Haver, 2014
		y=((m*wid/2+x-pos).^2)./(((wid/2).^2)+(x-pos).^2)
		% y=((1+(x-pos./(m.*wid))).^2)./(1+((x-pos)./wid).^2)
		g=y./max(y)
		% ----------------------------------------------------------------------
		def err = fitnbinpdf(tau,x,y)
		% Fitting def for iterative fit to the sum of
		% Negative Binomial Distributions
		% (http://www.mathworks.com/help/stats/negative-binomial-distribution.html)
		global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
		A = zeros(length(x),round(length(tau)/2))
		for j = 1:length(tau)/2,
			A(:,j) = nbinpdf(x,tau(2*j-1),tau(2*j))

		if AUTOZERO==3,A=[ones(size(y))' A]
		if BIPOLAR,PEAKHEIGHTS=A\y'else: PEAKHEIGHTS=abs(A\y')
		z = A*PEAKHEIGHTS
		if LOGPLOT,
			err = norm(log10(z)-log10(y)')
		else:
			err = norm(z-y')

		% ----------------------------------------------------------------------
		def err = fitlognpdf(tau,x,y)
		% Fitting def for iterative fit to the sum of
		% Lognormal Distributions
		% (http://www.mathworks.com/help/stats/lognormal-distribution.html)
		global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
		A = zeros(length(x),round(length(tau)/2))
		for j = 1:length(tau)/2,
			A(:,j) = lognormal(x,tau(2*j-1),tau(2*j))

		if AUTOZERO==3,A=[ones(size(y))' A]
		if BIPOLAR,PEAKHEIGHTS=A\y'else: PEAKHEIGHTS=abs(A\y')
		z = A*PEAKHEIGHTS
		if LOGPLOT,
			err = norm(log10(z)-log10(y)')
		else:
			err = norm(z-y')

		% ----------------------------------------------------------------------
		def g = lognormal(x,pos,wid)
		% lognormal def.  pos=position wid=half-width (both scalar)
		% lognormal(x,pos,wid), where x may be scalar, vector, or matrix
		% pos=position wid=half-width (both scalar)
		% T. C. O'Haver, 1991  
		g = exp(-(log(x/pos)/(0.01.*wid)) .^2)
		% ----------------------------------------------------------------------
		def err = fitsine(tau,x,y)
		% Fitting def for iterative fit to the sum of
		% sine waves (alpha test, NRFPT)
		global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
		A = zeros(length(x),round(length(tau)/2))
		for j = 1:length(tau)/2,
			A(:,j) = sine(x,tau(2*j-1),tau(2*j))

		if AUTOZERO==3,A=[ones(size(y))' A]
		if BIPOLAR,PEAKHEIGHTS=A\y'else: PEAKHEIGHTS=abs(A\y')
		z = A*PEAKHEIGHTS
		if LOGPLOT,
			err = norm(log10(z)-log10(y)')
		else:
			err = norm(z-y')

		% ----------------------------------------------------------------------
		def g=sine(x,f,phase) 
		% Sine wave (alpha test)
		g=sin(2*pi*f*(x+phase))
		% ----------------------------------------------------------------------
		def err = fitd1gauss(lambda,t,y)
		%   Fitting defs for the first derivative of a Gaussian
		%  T. C. O'Haver, 2014
		global PEAKHEIGHTS AUTOZERO BIPOLAR
		A = zeros(length(t),round(length(lambda)/2))
		for j = 1:length(lambda)/2,
			A(:,j) = d1gauss(t,lambda(2*j-1),lambda(2*j))'

		if AUTOZERO==3,A=[ones(size(y))' A]
		if BIPOLAR,PEAKHEIGHTS=A\y'else: PEAKHEIGHTS=abs(A\y')
		z = A*PEAKHEIGHTS
		err = norm(z-y')
		% ----------------------------------------------------------------------
		def y=d1gauss(x,p,w)
		% First derivative of Gaussian (alpha test)
		y=-(5.54518.*(x-p).*exp(-(2.77259.*(p-x).^2)./w^2))./w.^2
		y=y./max(y)
		% ----------------------------------------------------------------------
		def coeff = fitpolynomial(t,y,order)
		coeff=polyfit(t,y,order)
		% order=order
		% coeff=coeff
		% ----------------------------------------------------------------------
		def y=polynomial(t,coeff)
		y=polyval(coeff,t)
		% ----------------------------------------------------------------------
		def err = fitsegmented(lambda,t,y)
		%   Fitting defs for articulated segmented linear
		%  T. C. O'Haver, 2014
		global LOGPLOT
		breakpoints=[t(1) lambda max(t)]
		z = segmented(t,y,breakpoints)
		% lengthz=length(z)
		if LOGPLOT,
			err = norm(log10(z)-log10(y)')
		else:
			err = norm(z-y)

		% ----------------------------------------------------------------------
		def yi=segmented(x,y,segs)
		global PEAKHEIGHTS
		clear yy
		for n=1:length(segs)
		yind=val2ind(x,segs(n))
		yy(n)=y(yind(1))

		yi=INTERP1(segs,yy,x)
		PEAKHEIGHTS=segs
		% ----------------------------------------------------------------------
		def err = fitlinslope(tau,x,y)
		% Fitting def for iterative fit to linear def
		global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
		A = zeros(length(x),round(length(tau)/2))
		for j = 1:length(tau)/2,
			z = (x.*tau(2*j-1)+tau(2*j))'
			A(:,j) = z./max(z)

		if AUTOZERO==3,A=[ones(size(y))' A]
		if BIPOLAR,PEAKHEIGHTS=A\y'else: PEAKHEIGHTS=abs(A\y')
		z = A*PEAKHEIGHTS
		if LOGPLOT,
			err = norm(log10(z)-log10(y)')
		else:
			err = norm(z-y')

		% ----------------------------------------------------------------------
		def y=linslope(x,slope,intercept)
		y=x.*slope+intercept
		% y=y./max(y)
		% ----------------------------------------------------------------------
		def b=iqr(a)
		% b = IQR(a)  returns the interquartile range of the values in a.  For
		%  vector input, b is the difference between the 75th and 25th percentiles
		%  of a.  For matrix input, b is a row vector containing the interquartile
		%  range of each column of a.
		%  T. C. O'Haver, 2012
		mina=min(a)
		sizea=size(a)
		NumCols=sizea(2)
		for n=1:NumCols,b(:,n)=a(:,n)-mina(n)
		Sorteda=sort(b)
		lx=length(Sorteda)
		SecondQuartile=round(lx/4)
		FourthQuartile=3*round(lx/4)
		b=abs(Sorteda(FourthQuartile,:)-Sorteda(SecondQuartile,:))
		% ----------------------------------------------------------------------
		def err = fitmultiple(lambda,t,y,shapesvector,m)
		% Fitting def for a multiple-shape band signal.
		% The sequence of peak shapes are defined by the vector "shape".
		% The vector "m" determines the shape of variable-shape peaks.
		global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT coeff
		numpeaks=round(length(lambda)/2)
		A = zeros(length(t),numpeaks)
		for j = 1:numpeaks,
			if shapesvector(j)==28,
					coeff=polyfit(t,y,m(j))
					A(:,j) = polyval(coeff,t)
			else:
					A(:,j) = peakdef(shapesvector(j),t,lambda(2*j-1),lambda(2*j),m(j))'
			

		if AUTOZERO==3,A=[ones(size(y))' A]
		if BIPOLAR,PEAKHEIGHTS=A\y'else: PEAKHEIGHTS=abs(A\y')
		z = A*PEAKHEIGHTS
		if LOGPLOT,
			err = norm(log10(z)-log10(y)')
		else:
			err = norm(z-y')

		% ----------------------------------------------------------------------
		def p=peakdef(shape,x,pos,wid,m,coeff)
		% def that generates any of 20 peak types specified by number. 'shape'
		% specifies the shape type of each peak in the signal: "peakshape" = 1-20.
		% 1=Gaussian 2=Lorentzian, 3=logistic, 4=Pearson, 5=exponentionally
		% broadened Gaussian 9=exponential pulse, 10=up sigmoid,
		% 13=Gaussian/Lorentzian bl 14=BiGaussian, 15=Breit-Wigner-Fano (BWF) ,
		% 18=exponentionally broadened Lorentzian 19=alpha def 20=Voigt
		% profile 21=triangular 23=down sigmoid 25=lognormal. "m" is required
		% for variable-shape peaks only.
		switch shape,
			case 1
					p=gaussian(x,pos,wid)
			case 2
					p=lorentzian(x,pos,wid)
			case 3
					p=logistic(x,pos,wid)
			case 4
					p=pearson(x,pos,wid,m)
			case 5
					p=expgaussian(x,pos,wid,m)
			case 6
					p=gaussian(x,pos,wid)
			case 7
					p=lorentzian(x,pos,wid)
			case 8
					p=expgaussian(x,pos,wid,m)'
			case 9
					p=exppulse(x,pos,wid)
			case 10
					p=upsigmoid(x,pos,wid)
			case 11
					p=gaussian(x,pos,wid)
			case 12
					p=lorentzian(x,pos,wid)
			case 13
					p=GL(x,pos,wid,m)
			case 14
					p=BiGaussian(x,pos,wid,m)
			case 15
					p=BWF(x,pos,wid,m)
			case 16
					p=gaussian(x,pos,wid)
			case 17
					p=lorentzian(x,pos,wid)
			case 18
					p=explorentzian(x,pos,wid,m)'
			case 19
					p=alphadef(x,pos,wid)
			case 20
					p=voigt(x,pos,wid,m)
			case 21
					p=triangular(x,pos,wid)    
			case 23
					p=downsigmoid(x,pos,wid)
			case 25
					p=lognormal(x,pos,wid)
			case 26
					p=linslope(x,pos,wid)
			case 27
					p=d1gauss(x,pos,wid)
			case 28
					p=polynomial(x,coeff)
			otherwise
		% switch
