# Exercise 2.3-b - Error Euler-Cromer method

import math, pylab
g = 9.8 # acceleration due to gravity (m/s^2)
l = input("Length of pendulum (meters) ->")
theta0 = input("Initial angle (radians) ->")
dtlow = input("Lowest resolution time step (seconds)->")
nres = input("Number of resolution refinements ->")
tmax = input("Time to end of simulation (seconds) ->")
omega0 = math.sqrt(g/l)
pylab.figure()
for n in range(nres):
    refine = 10**n
    dt = dtlow/refine
    nsteps = int(tmax/dt)
    t = [0.0]*nsteps
    # omega, theta and error for Euler-Cromer method
    omega = 0.0
    theta = theta0
    err = [0.0] * nsteps
    for i in range(nsteps -1):
        t[i+1] = t[i]+dt
        #  Euler-Cromer method
        omega = omega -(g/l)*theta*dt
        theta = theta + omega * dt
        # Errors
        exact = theta0*math.cos(omega0*t[i+1])
        err[i+1] = abs(theta-exact)
    # Plotting the errors at this resolution
    pylab.loglog(t[refine::refine],
    err[refine::refine],'.-',label='dt = ' +str(dt))
pylab.legend(loc=3)
pylab.xlabel('Time - log scale')
pylab.ylabel('Absolute Error - log scale')
pylab.title('Euler-Cromer method integration error')
pylab.grid(linestyle='-', which='major')
pylab.grid(which='minor')
pylab.show()







