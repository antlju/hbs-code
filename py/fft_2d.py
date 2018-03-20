from numpy import *
from matplotlib import pyplot as plt

N = 4
L0 = -2.0*pi
L1 = 2.0*pi
x = linspace(L0,L1,N+1)
y = x
q1 = 1.0
q2 = 1.0

A = zeros((len(x),len(x)))
Aa = zeros((len(x),len(x)))
dimlen = arange(0,len(x))
for i in dimlen:
    for j in dimlen:
        A[i][j] = cos(q1*x[i])*cos(q2*y[j])
        Aa[i][j] = -cos(q1*x[i])*cos(q2*y[j])/(q1**2+q2**2)


Ahat = fft.rfft2(A)
Aback = fft.irfft2(Ahat,A.shape)

dx = x[1]-x[0]
for i in arange(0,Ahat.shape[0]):
    for j in arange(0,Ahat.shape[1]):
        fac = power(2*pi*(i*dx),2)+power(2*pi*(j*dx),2)
        if fac == 0:
            Ahat[i][j] = 0
        else:
            Ahat[i][j] = (-1.0/fac)*Ahat[i][j]
            

Psi = fft.irfft2(Ahat,A.shape)

#print absolute(Psi-Aa)
print absolute(Aback-A)


