# Exercise 2.5-d - SHM Euler-Cromer - Error at different resolutions

import math, pylab
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
    # Euler-Cromer method
    for i in range(nsteps - 1):
        v = v - B*x*dt
        x = x + v*dt
        t[i+1]= t[i] + dt

        exact = v0*math.sin(math.sqrt(B)*t[i+1])+x0*math.cos(math.sqrt(B)*t[i+1])
        err[i+1] = abs(x-exact)
    pylab.loglog(t[refine::refine],err[refine::refine],'.-',label='dt = '+str(dt))
pylab.legend(loc = 2)
pylab.xlabel('Time - log scale')
pylab.ylabel('Absolute Error - log scale')
pylab.title('Euler-Cromer method integration error')
pylab.grid(linestyle='-', which='major')
pylab.grid(which='minor')
pylab.show()

