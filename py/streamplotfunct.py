from numpy import *
from matplotlib import pyplot as plt


### This script plots the x-component of u(x,y,z) in the (x,y)-plane.

### Load data. This is a 1D array using C-style access with last index being the fastest.
#fname0 = "../simdata/kolm_rk3_x_N64_tNum_60.dat"
def kolmstreamplot(fn,nn,kf):
    fname0 = fn   

    u0 = loadtxt(fname0,skiprows=1)
    N = nn 
    Nx = N
    Ny = N
    Nz = N
    
    U0 = zeros((N,N,N));
    
    for i in range(0,N):
        for j in range(0,N):
            for k in range(0,N):
                U0[i][j][k] = u0[(Ny*i+j)*Nz+k]
                
                # Fix k=0

    #print U0[:-2]
    Umesh = U0[:][:][0]

    #print(Umesh)

    x = linspace(-pi,pi,N)
    y = x

    X,Y = meshgrid(x,y)

    Vzer = zeros((N,N))
    #skip=(slice(None,None,4),slice(None,None,4))
    plt.figure(1)
    #plt.quiver(x[skip],y[skip],Umesh[skip],Vzer[skip],Umesh[skip])
    plt.quiver(x,y,Umesh,Vzer,Umesh)
    #plt.imshow(Umesh)
    tit = "kf = %i, N = %i" % (kf,nn)
    plt.title(tit)
    plt.xlabel('x')
    plt.ylabel('y')
    
    plt.show()

def meanEplot(fname):
    #meanE = loadtxt(fname,usecols=[3])
    meanE = loadtxt(fname)
    #timesteps = loadtxt(fname,usecols=0,delimiter=":")

    plt.plot(meanE)
    plt.show()
