# Phy-801
# Exercise 3.1-a - leap-frog method

import math
import pylab
#fixed parameters
c = 300.0 #speed of propagation (m/s)
L = 1.0 #length of wire
x0 = 0.4 #initial pulse location
s = 0.02 #initial pulse width
#input parameters
##dx = L * input('Grid spacing in units of wire length (L) ->')
##dt = (dx/c) * input('Time step in units of (dx/c) ->')
##tmax = (L/c) * input('Evolution time in units of (L/c) ->')
dx = 0.01        ##(for input above 0.01)
dt = (dx/c)  ##(for input above 1.0, 0.5)
tmax = 0.0025     ##(for input above 1.2)
#construct initial data
N = int(L / dx)
x = [0.0] * (N+1)
u0 = [0.0] * (N+1)
v0 = [0.0] * (N+1)
u1 = [0.0] * (N+1)
v1 = [0.0] * (N+1)
for j in range(N+1):
    x[j] = j * dx
    u0[j] = math.exp(-0.5*((x[j]-x0)/s)**2)
#perform the evolution
t = 0.0
p = 0
tstep = 0.0003
while t < tmax:
    #derivatives at interior points (leap-frog method)
    for j in range(N):
        v1[j] = v0[j] + dt * c * (u0[j+1] - u0[j])/dx   
    for j in range(1,N):
        u1[j] = u0[j] + dt * c * (v1[j] - v1[j-1])/dx
    #boundary conditions
    u1[0] = u1[1]
    u1[N] = u1[N-1]    
    #plotting
    if t >= p*tstep:
                pylab.plot(x, [u+1.1*p+1  for u in u0], label='t (s)='+str(round(t,4)) )
                pylab.legend(loc=(1.03,0.2))
                p += 1
    #swap old and new lists
    u0, u1 = u1, u0
    v0, v1 = v1, v0
    t += dt
pylab.ylim(0.0,11.0)
pylab.xlabel('x(m)')
pylab.ylabel('u')
pylab.title('Evolution of wave using leap-frog method')
pylab.show()



















    

