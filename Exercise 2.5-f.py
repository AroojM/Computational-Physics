# Exercise 2.5-f - SHM RK4 - Error at different resolutions

import math, pylab
# Defining functions that return values of velocity and acceleration
def f1(t, x, v):
    return v
def f2(t, x, v):
    a = -B*x
    return a
x0 = input('Initial x position (m) ->')
v0 = input('Initial velocity (m/s)->')
B= input('Constant in harmonic motion equation (1/s) ->')
dtlow = input('Lowest resolution time step (s) ->')
nres = input('Number of resolution refinements ->')
tmax = input('Time to the end of integration (s) ->')
pylab.figure()
for n in range(nres):
    refine = 10**n
    dt = dtlow/refine
    nsteps = int(tmax / dt)
    t = [0.0] * nsteps
    err = [0.0] * nsteps
    x = x0
    v = v0
    # rk4 method
    for i in range(nsteps - 1):
        K1_1=f1(t[i],x,v)
        K1_2=f2(t[i],x,v)
        K2_1=f1(t[i]+0.5*dt,x+0.5*dt*K1_1,v+0.5*dt*K1_2)
        K2_2=f2(t[i]+0.5*dt,x+0.5*dt*K1_1,v+0.5*dt*K1_2)
        K3_1=f1(t[i]+0.5*dt,x+0.5*dt*K2_1,v+0.5*dt*K2_2)
        K3_2=f2(t[i]+0.5*dt,x+0.5*dt*K2_1,v+0.5*dt*K2_2)
        K4_1=f1(t[i]+dt,x+dt*K3_1,v+dt*K3_2)
        K4_2=f2(t[i]+dt,x+dt*K3_1,v+dt*K3_2)

        x = x + (dt/6.0)*(K1_1+2.0*K2_1+2.0*K3_1+K4_1)
        v = v + (dt/6.0)*(K1_2+2.0*K2_2+2.0*K3_2+K4_2)
        t[i+1]= t[i] + dt

        exact = v0*math.sin(math.sqrt(B)*t[i+1])+x0*math.cos(math.sqrt(B)*t[i+1])
        err[i+1] = abs(x-exact)
    pylab.loglog(t[refine::refine],err[refine::refine],'.-',label='dt = '+str(dt))
pylab.legend(loc = 2)
pylab.xlabel('Time - log scale')
pylab.ylabel('Absolute Error - log scale')
pylab.title('RK4 method integration error')
pylab.grid(linestyle='-', which='major')
pylab.grid(which='minor')
pylab.show()

































