from numpy import *
from matplotlib import pyplot as plt
import matplotlib.animation as animation

fname0 = "/home/anton/dev/hbs/simdata/kolmo_v1_component_0__tNum_99.csv"
fname1 = "/home/anton/dev/hbs/simdata/kolmo_v1_component_1__tNum_99.csv"
fname2 = "/home/anton/dev/hbs/simdata/kolmo_v1_component_2__tNum_99.csv"

F0 = loadtxt(fname0,skiprows=1)
F1 = loadtxt(fname1,skiprows=1)  



N = 32
F0 = F0[0:N]
F1 = F1[0:N]

L0 = 0
L1 = 2*pi
x = linspace(L0,L1-L1/N,N)
y = x

X,Y = meshgrid(x,y)

###plt.streamplot(X,Y,F0,F1,color=F1, linewidth=2, cmap='autumn',density=[0.5,1])
#plt.quiver(F0,F1)
plt.imshow(F0)
plt.show()
