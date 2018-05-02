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

fname = "/home/anton/dev/hbs/simdata/rk4test_gaussInit_2_tNum_0.csv"

X = loadtxt(fname, usecols=3,delimiter=",",skiprows=1)
Y = loadtxt(fname, usecols=4,delimiter=",",skiprows=1)
Z = loadtxt(fname, usecols=5,delimiter=",",skiprows=1)
F = loadtxt(fname, usecols=6,delimiter=",",skiprows=1)


R2 = power(X,2)+power(Y,2)+power(Z,2)
R2,Ridx = unique(R2,return_index=True)
R2 = sort(R2)
R = sqrt(R2)
F = F[Ridx]

lenF = len(F)
lenR = len(R)

Fbig = F
Rbig = R

def gaussian(rr):
    prefac = 1#/power(2*pi,3)
    ret = prefac*exp(-power(rr,2))

    return ret

ims = []

fig = plt.figure()
ax = plt.axes(xlim=(0, 10), ylim=(0, 1))
dt = 0.0061685
#num = 99
for num in arange(0,200,10):
    fname = "/home/anton/dev/hbs/simdata/rk4test_deltaInit_2_tNum_"+str(num)+".csv"

    X = loadtxt(fname, usecols=3,delimiter=",",skiprows=1)
    Y = loadtxt(fname, usecols=4,delimiter=",",skiprows=1)
    Z = loadtxt(fname, usecols=5,delimiter=",",skiprows=1)
    F = loadtxt(fname, usecols=6,delimiter=",",skiprows=1)

    N = len(X)
    R2 = power(X,2)+power(Y,2)+power(Z,2)
    #R2 = R2[len(R2)/2:-1]
    Nr = int(sqrt(3)*N/2)
    den_state = zeros(Nr+1)
    Ur = zeros(Nr+1)
    #Rr = zeros(Nr+1)

    #print R2
    for idx in arange(0,len(R2)):
        ir = int(sqrt(R2[idx]))
        Ur[ir] = Ur[ir] + F[idx]
        #Rr[ir] = ir
        den_state[ir] = den_state[ir] + 1

    for i in arange(len(Ur)):
        if den_state[i] == 0:
            Ur[i] = 0
        else:
            Ur[i] = Ur[i]/den_state[i]
            
    ims.append(plt.plot(Ur/sum(Ur),'b.-',linewidth=2))
    #plt.plot(Ur/sum(Ur),'b.-',linewidth=2)
    #Rr = linspace(0,10,10)
    #plt.plot(u_d(Rr,num*dt),'r.-',linewidth=2)
        

    

    

im_ani = animation.ArtistAnimation(fig, ims, interval=150, repeat_delay=200,
                                   blit=True)

plt.show()
