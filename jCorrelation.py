#jCorrelation.py
from numpy.fft import fft,ifft
from numpy import loadtxt,savetxt
a=loadtxt("jin.txt");
b=fft(a[:,1])
b1=b.real;
b2=b.imag;
c=b1*b1+b2*b2;
d=ifft(c);
savetxt("jcor.txt",d);