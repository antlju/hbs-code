from numpy import *
from matplotlib import pyplot as plt
import matplotlib.animation as animation

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
ax = plt.axes(xlim=(0, 15), ylim=(0, 1))

#num = 99
for num in arange(0,200,1):
    fname = "/home/anton/dev/hbs/simdata/rk4test_deltaInit_2_tNum_"+str(num)+".csv"

    X = loadtxt(fname, usecols=3,delimiter=",",skiprows=1)
    Y = loadtxt(fname, usecols=4,delimiter=",",skiprows=1)
    Z = loadtxt(fname, usecols=5,delimiter=",",skiprows=1)
    F = loadtxt(fname, usecols=6,delimiter=",",skiprows=1)


    R2 = power(X,2)+power(Y,2)+power(Z,2)
    R2,Ridx = unique(R2,return_index=True)
    R2 = sort(R2)
    R = sqrt(R2)
    F = F[Ridx]
    #Fbig = vstack((Fbig,F))
    ims.append(plt.plot(R,F,'b.-',linewidth=2))

im_ani = animation.ArtistAnimation(fig, ims, interval=150, repeat_delay=200,
                                   blit=True)

plt.show()

#plt.plot(R,Fbig[0])
#plt.show()











#fig = plt.figure()
#ax = plt.axes(xlim=(0, 10), ylim=(0, 1))
#line, = ax.plot([], [], lw=2)

#def init():
#    line.set_data([], [])
#    return line,

#def animate(i,Rbig,Fbig):
#    x = Rbig[i]
#    y = Fbig[i]
#    line.set_data(x, y)
#    return line,

#anim = animation.FuncAnimation(fig, animate, frames=len(Fbig), fargs=(Rbig,Fbig),interval=100, blit=False)

#anim.save('basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

#plt.show()
