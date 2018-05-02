from numpy import *
from matplotlib import pyplot as plt
import matplotlib.animation as animation

fname0 = "/home/anton/dev/hbs/simdata/kolmo_v1_component_0__tNum_49.csv"
fname1 = "/home/anton/dev/hbs/simdata/kolmo_v1_component_1__tNum_49.csv"
fname2 = "/home/anton/dev/hbs/simdata/kolmo_v1_component_2__tNum_49.csv"

X = loadtxt(fname0, usecols=3,delimiter=",",skiprows=1)
Y = loadtxt(fname0, usecols=4,delimiter=",",skiprows=1)
Z = loadtxt(fname0, usecols=5,delimiter=",",skiprows=1)
F0 = loadtxt(fname0, usecols=6,delimiter=",",skiprows=1)
F1 = loadtxt(fname1, usecols=6,delimiter=",",skiprows=1)
F2 = loadtxt(fname2, usecols=6,delimiter=",",skiprows=1)   

N = 8

L0 = 0
L1 = 2*pi
x = linspace(L0,L1-L1/N,N)
y = x

MX,MY = meshgrid(X,Y)

for j in arange(0,N):
    for k in arange(0,N):
        print F0[k+j*N]

MF0 = meshgrid(F0[0:(N*N):N])

print MF0
