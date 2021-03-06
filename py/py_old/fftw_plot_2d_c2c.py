from numpy import *
from matplotlib import pyplot as plt
import matplotlib.animation as animation


fname1 = "/home/anton/dev/hbs/src/cosxcosy_2d.out"
fname3 = "/home/anton/dev/hbs/src/sinxsiny_2d.out"
fname2 = "/home/anton/dev/hbs/src/4diracs.out"

N=64

FFT = loadtxt(fname3,delimiter=',',usecols=range(N))

L0 = 0
L1 = 2*pi
x = linspace(L0,L1,N)
y = x

X,Y = meshgrid(x,y)

q1 = 1.0
q2 = q1
#analytic = -sin(q1*X)*sin(q2*Y)/(power(q1,2)+power(q2,2))
analytic = zeros((N,N))
icanalytic = zeros((N,N))
for i in range(N):
    for j in range(N):
        #analytic[i][j] = -cos(q1*x[i])*cos(q2*y[j])/(power(q1,2)+power(q2,2))
        #icanalytic[i][j] = cos(q1*x[i])*cos(q2*y[j])
        analytic[i][j] = -sin(q1*x[i])*sin(q2*y[j])/(power(q1,2)+power(q2,2))

print absolute(analytic-FFT)

plt.figure(1)
plt.imshow(FFT)

plt.figure(2)
plt.imshow(analytic)
plt.show()
