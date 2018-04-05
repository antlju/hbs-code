from numpy import *
from matplotlib import pyplot as plt
import matplotlib.animation as animation



fname1 = "/home/anton/dev/hbs/src/sinxtest.out"

FFT = loadtxt(fname1)

N = len(FFT)

L0 = 0
L1 = 2*pi
x = linspace(L0,L1-L1/N,N)
print x

A = sin(x)
plt.plot(x,FFT)
#plt.plot(x,A)
plt.show()
