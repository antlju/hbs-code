from numpy import *
from matplotlib import pyplot as plt


### This script plots the x-component of u(x,y,z) in the (x,y)-plane.

### Load data. This is a 1D array using C-style access with last index being the fastest.
#fname0 = "../simdata/kolm_rk3_x_N64_tNum_60.dat"
#fname0 = "../simdata/kolmo_kf_1_N_256_stepNo_40.dat"   
fname0 = "../simdata/kolmo_kf_1_N_128_stepNo_61.dat"

u0 = loadtxt(fname0,skiprows=1)
#N = 256
N = 128
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

#print Umesh

x = linspace(-pi,pi,N)
y = x

X,Y = meshgrid(x,y)

Vzer = zeros((N,N))
plt.figure(1)
plt.quiver(x,y,Umesh,Vzer,Umesh+0.5)
#plt.imshow(Umesh)
plt.xlabel('x')
plt.ylabel('y')
    
plt.show()


