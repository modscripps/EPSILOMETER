function kaxis=make_kaxis(tscan,df)

% create the frequency axis for fft
% tscan is the length in time units of the segment for the fft
% df is the sampling frequency in time units^{-1}
Lscan  = floor(tscan*df);
dk=1/tscan;
if rem(Lscan,2)==0
    kaxis=-Lscan/2*dk:dk:Lscan/2*dk-dk;
else
    kp=dk:dk:dk*floor(Lscan/2);
    kaxis=[fliplr(-kp) 0 kp];
end
kaxis=fftshift(kaxis);


