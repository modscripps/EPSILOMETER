function z2=remove_sbe_preemphasis(data,fs);
%function z2=remove_sbe_preemphasis(data,fs);
%Remove preemphasis filter, following Steve Anderson and Alfor/Pinkel 2000,
%from SBE conductivity sensor.

%This is taken from Steve's preemphasis2.
i=sqrt(-1);
clear G;
clear X  y  Y;
PI=3.14159;
R24=1e6;
R25=577e3;
Rf=R25+R24;
R22=266.1;
C19=1e-6;
if mod(length(data),2)~=0
    n=length(data(1:end-1));
    data=data(1:end-1);
else
    n=length(data);
end
%w=0:96/n:(48-96/n);
%fs=160; %was 96
w=0:fs/n:(fs/2-fs/n);

h=1;
w1=h./(w*PI*2*C19);
G=(1+Rf*R22./(R22*R22+w1.^2)) + i*(Rf*w1./(R22*R22+w1.^2));

G=h./G;
G(1)=sqrt(2);

n=length(G);
m=2*n;
G(n+1)=abs(G(n));
G(n+2:m)=conj(G(n:-1:2));
G(1)=0;
%z=FCTD.ucon;


y=fft(detrend(data));
y=y.*G.';

z2=real(ifft(y));
