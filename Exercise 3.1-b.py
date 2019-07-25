# Exercise 3.1-b - wave-lax method

import math
import pylab
#fixed parameters
c = 300.0 #speed of propagation (m/s)
L = 1.0 #length of wire
x0 = 0.3 #initial pulse location
s = 0.02 #initial pulse width
#input parameters
## A = input('Amplitude of sinusoid ->')
## omega = input('Frequency of sinusoid ->')
##dx = L * input('Grid spacing in units of wire length (L) ->')
##dt = (dx/c) * input('Time step in units of (dx/c) ->')
##tmax = (L/c) * input('Evolution time in units of (L/c) ->')
# input values: dx = 0.01, dt = 1.0, tmax = 15
dx = 0.01
dt =  (dx/c)
tmax = 0.05
A = 1.0
omega = (10.0*math.pi*c)/(L)
#construct initial data
N = int(L / dx)
x = [0.0] * (N+1)
u0 = [0.0] * (N+1)
v0 = [0.0] * (N+1)
u1 = [0.0] * (N+1)
v1 = [0.0] * (N+1)
for j in range(N+1):
    x[j] = j * dx
#prepare animated plot
pylab.ion()
line, = pylab.plot(x, u0, '-k')
pylab.ylim(-4.0,4.0)
pylab.xlabel('x(m)')
pylab.ylabel('u')
#perform the evolution
t = 0.0
while t < tmax:
    #update plot
    line.set_ydata(u0)
    pylab.title('t = %5f' % t)
    pylab.draw()
    pylab.pause(0.001)
    #derivatives at interior points
    for j in range(1,N):
        v1[j] = 0.5*(v0[j-1]+v0[j+1]) + 0.5 * dt * c * (u0[j+1] - u0[j-1])/dx
        u1[j] = 0.5*(u0[j-1]+u0[j+1]) + 0.5 * dt * c * (v0[j+1] - v0[j-1])/dx       
    #boundary conditions
    u1[0] =A*math.sin(omega*t)
    u1[N] = 0.0
    v1[0] = v1[1] + u1[1]
    v1[N] = v1[N-1] - u1[N-1]
    #swap old and new lists
    u0, u1 = u1, u0
    v0, v1 = v1, v0
    t += dt
#freeze final plot
pylab.ioff()
pylab.show()




















    

