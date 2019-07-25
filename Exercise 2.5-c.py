# Exercise 2.5-c - SHM RK4
# Exercise 2.5-c

import math
import pylab
# Defining functions that return values of velocity and acceleration
def f1(t, x, v):
    return v
def f2(t, x, v):
    a = -B*x
    return a
x0 = input('Initial x position (m) ->')
v0 = input('Initial velocity (m/s)->')
dt = input('Time step (s) ->')
tmax = input('Time to the end of integration (s) ->')
B= input('Constant in harmonic motion equation (1/s) ->')
nsteps = int(tmax / dt)
x = [0.0] * nsteps
t = [0.0] * nsteps
t[0] = 0
x[0] = x0
v = v0
# Integrating by rk4 method
for i in range(nsteps - 1):
    K1_1=f1(t[i],x[i],v)
    K1_2=f2(t[i],x[i],v)
    K2_1=f1(t[i]+0.5*dt,x[i]+0.5*dt*K1_1,v+0.5*dt*K1_2)
    K2_2=f2(t[i]+0.5*dt,x[i]+0.5*dt*K1_1,v+0.5*dt*K1_2)
    K3_1=f1(t[i]+0.5*dt,x[i]+0.5*dt*K2_1,v+0.5*dt*K2_2)
    K3_2=f2(t[i]+0.5*dt,x[i]+0.5*dt*K2_1,v+0.5*dt*K2_2)
    K4_1=f1(t[i]+dt,x[i]+dt*K3_1,v+dt*K3_2)
    K4_2=f2(t[i]+dt,x[i]+dt*K3_1,v+dt*K3_2)
   
    x[i+1] = x[i] + (dt/6.0)*(K1_1+2.0*K2_1+2.0*K3_1+K4_1)
    v = v + (dt/6.0)*(K1_2+2.0*K2_2+2.0*K3_2+K4_2)
    t[i+1]= t[i] + dt
pylab.plot(t, x)
pylab.xlabel('t (s)')
pylab.ylabel('x (m)')
pylab.title('Simple harmonic oscillator - x VS t - RK4 method')
pylab.grid()
pylab.show()






























