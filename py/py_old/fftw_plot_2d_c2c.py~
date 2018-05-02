from numpy import *
from matplotlib import pyplot as plt
import matplotlib.animation as animation

fname = "/home/anton/dev/hbs/src/sinxsiny.out"
X = loadtxt(fname, usecols=0,delimiter="\t",skiprows=0)
Y = loadtxt(fname, usecols=1,delimiter="\t",skiprows=0)
F = loadtxt(fname, usecols=2,delimiter="\t",skiprows=0)

N = 16
FV = zeros((N,N))
for i in arange(0,N):
    for j in arange(0,N):
        FV[i][j] = F[j+N*i]
x = linspace(0,2*pi,N)
y = linspace(0,2*pi,N)
XV, YV = meshgrid(x,y)
Z = -cos(XV)*cos(YV)/2
plt.contourf(XV,YV,FV)
plt.show()
