from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from numpy import *
from matplotlib import pyplot as plt
import matplotlib.pylab as pl
from matplotlib.colors import ListedColormap

plt.figure()

num = 0
for num in arange(0,200,20):
    fname = "/home/anton/dev/hbs/simdata/rk2test_gaussInit_2_tNum_"+str(num)+".csv"

    I = loadtxt(fname, usecols=0,delimiter=",",skiprows=1)
    J = loadtxt(fname, usecols=1,delimiter=",",skiprows=1)
    K = loadtxt(fname, usecols=2,delimiter=",",skiprows=1)
    X = loadtxt(fname, usecols=3,delimiter=",",skiprows=1)
    Y = loadtxt(fname, usecols=4,delimiter=",",skiprows=1)
    Z = loadtxt(fname, usecols=5,delimiter=",",skiprows=1)
    F = loadtxt(fname, usecols=6,delimiter=",",skiprows=1)

    FF, FFidx = unique(F,return_index=True)
    R2 = pow(take(X,FFidx),2)+pow(take(Y,FFidx),2)+pow(take(Z,FFidx),2)
    R = sqrt(R2)
    R = sort(R)
    FF = flip(FF,0)

    
    plt.plot(R,FF,'-')
    
plt.show()
