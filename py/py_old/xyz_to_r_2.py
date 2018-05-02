from numpy import *
from matplotlib import pyplot as plt
import matplotlib.animation as animation

def u_d(R,t):
    if t == 0:
        ret = zeros(len(R))
        ret[0] = 1.0
    else:
        D = 1.0
        C = 1.0/(8*sqrt(power(pi,3)))
        C = 1.0
        ret = (C/(sqrt(power(t,3))))*exp(-power(R,2)/(4*D*t))

    return ret/sum(ret)

fname = "/home/anton/dev/hbs/simdata/rk4test_gaussInit_2_tNum_1.csv"

X = loadtxt(fname, usecols=3,delimiter=",",skiprows=1)
Y = loadtxt(fname, usecols=4,delimiter=",",skiprows=1)
Z = loadtxt(fname, usecols=5,delimiter=",",skiprows=1)
F = loadtxt(fname, usecols=6,delimiter=",",skiprows=1)

print len(X)
N = len(X)
R2 = power(X,2)+power(Y,2)+power(Z,2)
#R2 = R2[len(R2)/2:-1]
Nr = int(sqrt(3)*N/2)
den_state = zeros(Nr+1)
Ur = zeros(Nr+1)
Rr = zeros(Nr+1)

print R2
for idx in arange(0,len(R2)):
    ir = int(sqrt(R2[idx]))
    Ur[ir] = Ur[ir] + F[idx]
    den_state[ir] = den_state[ir] + 1
    

    
plt.figure()
plt.plot(Ur/den_state,'.-')
plt.show()
