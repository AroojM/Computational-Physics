# Exercise 2.5-a - SHM Euler-Cromer
# Exercise 2.5-a

import math
import pylab
x0 = input('Initial x position (m) ->')
v0 = input('Initial velocity (m/s)->')
dt = input('Time step (s) ->')
tmax = input('Time to the end of integration (s) ->')
B= input('Constant in harmonic motion equation (1/s) ->')
nsteps = int(tmax / dt)
x = [0.0] * nsteps
t = [0.0] * nsteps
# Initial conditions
t[0] = 0
x[0] = x0
v = v0
# Use Euler-Cromer method to integrate SNM equations
for i in range(nsteps - 1):
    v = v -B*x[i]*dt
    x[i+1] = x[i] + v*dt
    t[i+1] = t[i] + dt
pylab.plot(t, x)
pylab.xlabel('t (s)')
pylab.ylabel('x (m)')
pylab.title('Simple harmonic oscillator - x VS t - Euler-Cromer method')
pylab.grid()
pylab.show()






























