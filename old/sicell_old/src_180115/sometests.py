from numpy import *

x = linspace(0,3*pi/2,4)
y = x
z = y
k = 1.0

u1 = sin(k*z)+cos(k*y)
u2 = sin(k*x)+cos(k*z)
u3 = sin(k*y)+cos(k*x)

u = array([u1, u2, u3])

omega1 = k*cos(k*y)-(-k*sin(k*z))
omega2 = k*cos(k*z)-(-k*sin(k*x))
omega3 = k*cos(k*x)-(-k*sin(k*y))

omega = array([omega1, omega2, omega3])

udotu = u1*u1+u2*u2+u3*u3
udotw = u1*omega1+u2*omega2+u3*omega3

print udotw
print udotu
print udotw/udotu
