from numpy import *
from matplotlib import pyplot as plt
import matplotlib.animation as animation

def gaussian(rr):
    prefac = 1#/power(2*pi,3)
    ret = prefac*exp(-power(rr,2))

    return ret
    
num = 0
fname = "/home/anton/dev/hbs/simdata/rk4test_gaussInit_2_tNum_"+str(num)+".csv"

X = loadtxt(fname, usecols=3,delimiter=",",skiprows=1)
Y = loadtxt(fname, usecols=4,delimiter=",",skiprows=1)
Z = loadtxt(fname, usecols=5,delimiter=",",skiprows=1)
F = loadtxt(fname, usecols=6,delimiter=",",skiprows=1)


R2 = power(X,2)+power(Y,2)+power(Z,2)
R2,Ridx = unique(R2,return_index=True)
R2 = sort(R2)
R = sqrt(R2)
#print len(R), len(F[Ridx])

plt.plot(gaussian(R),'d-')
plt.plot(F[Ridx],'.-')
plt.xlim([0,30])
plt.show()
