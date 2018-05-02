from numpy import *
from matplotlib import pyplot as plt

N = 64
L0 = -2.0*pi
L1 = 2.0*pi
x = linspace(L0,L1,N)

q = 0.5
#a = zeros((N+1));
#a[N/2] = 1.0;
a = cos(q*x);

rho = fft.rfft(a)
#step = 2*L1/N
#k = 2*L1*fft.rfftfreq(N,d=step)
d = x[4]-x[3]
k = 2*pi*fft.rfftfreq(len(x), d)
k[0] = d
rho[0] = 0
Psihat = -rho/k
#print Psihat
Psi = fft.irfft(Psihat)
aa = -cos(q*x)/q#power(q,2)
plt.plot(x,Psi)
plt.plot(x,aa,'k--',linewidth=2)
plt.show()
