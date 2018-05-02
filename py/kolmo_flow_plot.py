from numpy import *
from matplotlib import pyplot as plt

### This script plots the x-component of u(x,y,z) in the (x,y)-plane.

### Load data. This is formatted as a set of (y,z)-arrays sequenced along x-direction.
fname0 = "/home/anton/dev/hbs/simdata/kolmo_v1_component_0_tNum_99.dat"
u0 = loadtxt(fname0,skiprows=1)
N = len(u0[0])


#print N

# The (x,y)-plane array
u0_xy = zeros((N,N))

### Set the (x,y)-plane array from loaded data.
for i in range(N):
    for j in range(N):
        u0_xy[i][j] = u0[N*i+j][0]

#print u0_xy


### Plot

plt.imshow(u0_xy.transpose()) #transpose so that we have a 'standard' y(x)-plot
plt.xlabel('x')
plt.ylabel('y')
plt.show()
