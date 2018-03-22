from numpy import *
from matplotlib import pyplot as plt

N = 4
L0 = 0.0#-2.0*pi
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
#Aback = fft.irfft2(Ahat,A.shape)

dimlen1 = arange(0,Ahat.shape[0])
dimlen2 = arange(0,Ahat.shape[1])
for i in dimlen1:
    for j in dimlen2:
        if i < len(x)/2:
            II = i
        else:
            II = len(x)-i

        fac = power(II,2)+power(j,2)
        if fac == 0:
            Ahat[i][j] = 0
        else:
            Ahat[i][j] = (1.0/fac)*Ahat[i][j]
        
Psi = fft.ihfft(Ahat)

print Psi-Aa
