# Exercise 3.1-c

import math
import pylab
#fixed parameters
c1 = 300.0 #speed of propagation (m/s)- left-hand side
c2 = 600.0 #speed of propagation (m/s)- right-hand side
L = 1.0 #length of wire
x0 = 0.4 #initial pulse location
s = 0.02 #initial pulse width
dx = 0.01
# choosing smaller dt based on speeds
if c1 > c2:
    dt =(dx/c1)
else:
    dt = (dx/c2)
tmax = 0.0009    
#construct initial data
N = int(L / dx)
x = [0.0] * (N+1)
u0 = [0.0] * (N+1)
v0 = [0.0] * (N+1)
u1 = [0.0] * (N+1)
v1 = [0.0] * (N+1)
for j in range(N+1):
    x[j] = j * dx
    u0[j] = 2.0*math.exp(-0.5*((x[j]-x0)/s)**2)
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
    pylab.title('Pulse in two strings: t (s) = %5f' % t)
    pylab.draw()
    pylab.pause(0.1)
    # leap frog method
    #derivatives at interior points
    for j in range(N):
        if j < N/2:
            v1[j] = v0[j] + dt * c1 * (u0[j+1] - u0[j])/dx
        else:
            v1[j] = v0[j] + dt * c2 * (u0[j+1] - u0[j])/dx
            
    for j in range(1,N):
        if j < N/2:
            u1[j] = u0[j] + dt * c1 * (v1[j] - v1[j-1])/dx
        else:
            u1[j] = u0[j] + dt * c2 * (v1[j] - v1[j-1])/dx
    #boundary conditions
    u1[0] = u1[N] = 0.0
    #swap old and new lists
    u0, u1 = u1, u0
    v0, v1 = v1, v0
    t += dt
#freeze final plot
pylab.ioff()
pylab.show()




















    

