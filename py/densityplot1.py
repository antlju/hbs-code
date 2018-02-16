from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from numpy import *
from matplotlib import pyplot as plt
import matplotlib.pylab as pl
from matplotlib.colors import ListedColormap

plt.figure()

#num = 5
for num in arange(0,200,20):
    #num = 100
    fname = "/home/anton/dev/hbs/src/simdata/u_tNum_"+str(num)+".dat"

    X = loadtxt(fname, usecols=0)
    Y = loadtxt(fname, usecols=1)
    Z = loadtxt(fname, usecols=2)
    U = loadtxt(fname, usecols=3)

    N = 32#len(X)
    
    Nr = int(sqrt(3)*N/2)+1
    Ur = zeros(Nr+1)
    den_state = zeros(Nr+1)
    #print Nr
    #print N
    for i in arange(0,N):
        for j in arange(0,N):
            for k in arange(0,N):
                r2 = (i-N/2)**2+(j-N/2)**2+(k-N/2)**2
                ir = int(rint(sqrt(r2)))#-sqrt(3)*16))
                #print ir
                Ur[ir] = Ur[ir] + U[i+N*(j+N*k)]
                den_state[ir] = den_state[ir] + 1
            
                #U[N*j+i(i+NG)+(N*j)+(N*k)
            
        
    for ir in arange(0,Nr):
        if den_state[ir] !=0:
            Ur[ir] = Ur[ir]/den_state[ir]
        else:
            Ur[ir] = 0
                        
                #print Ur
                
                
                    
    #print num
    
    plt.plot(Ur,'.-')
    plt.show()
    plt.hold(True)
    
    
                #print Ur



#plt.plot(Ur)

