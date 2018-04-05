from numpy import *
from matplotlib import pyplot as plt
import matplotlib.animation as animation

fname6 = "/home/anton/dev/hbs/src/sinxsiny_inplace.out"
FFT = loadtxt(fname6)
#FFT = FFT.transpose()
N = 64
L0 = 0
L1 = 2*pi
x = linspace(L0,L1-L1/N,N)
y = x

X,Y = meshgrid(x,y)

q1 = 2.0
q2 = 3.0
#analytic = -sin(q1*X)*sin(q2*Y)/(power(q1,2)+power(q2,2))
analytic = zeros((N,N))
icanalytic = zeros((N,N))
for i in range(N):
    for j in range(N):
        #analytic[i][j] = -cos(q1*x[i])*cos(q2*y[j])/(power(q1,2)+power(q2,2))
        #icanalytic[i][j] = cos(q1*x[i])*cos(q2*y[j])
        #analytic[i][j] = -sin(q1*x[i])*sin(q2*y[j])/(power(q1,2)+power(q2,2))
        analytic[i][j] = sin(q1*x[i])*sin(q2*y[j])/(power(q1,2)+power(q2,2))

plt.figure(1)
plt.imshow(FFT)

plt.figure(2)
plt.imshow(analytic)

plt.show()
