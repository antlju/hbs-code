from numpy import *
from matplotlib import pyplot as plt

### This script plots the x-component of u(x,y,z) in the (x,y)-plane.

### Load data. This is a 1D array using C-style access with last index being the fastest.
fname0 = "/home/anton/dev/hbs/simdata/kolm_rk3_x_tNum_30.dat"

u0 = loadtxt(fname0,skiprows=1)
N = 32
Nx = N
Ny = N
Nz = N

U0 = zeros((N,N,N));

for i in range(0,N):
    for j in range(0,N):
        for k in range(0,N):
            U0[i][j][k] = u0[(Ny*i+j)*Nz+k]
            
# Fix k=0

Umesh = U0[:][:][0]

plt.figure(1)
plt.imshow(Umesh)
plt.xlabel('x')
plt.ylabel('y')

plt.show()
