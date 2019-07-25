# Exercise 2.3-c - Error - Euler method

import math, pylab
g = 9.8 # acceleration due to gravity (m/s^2)
l = 1.0
theta0 = 0.2
dtlow = 0.1
nres = 4
tmax = 5
##l = input("Length of pendulum (meters) ->")
##theta0 = input("Initial angle (radians) ->")
##dtlow = input("Lowest resolution time step (seconds)->")
##nres = input("Number of resolution refinements ->")
##tmax = input("Time to end of simulation (seconds) ->")
omega0 = math.sqrt(g/l)
pylab.figure()
for n in range(nres):
    refine = 10**n
    dt = dtlow/refine
    nsteps = int(tmax/dt)
    t = [0.0]*nsteps
    # omega, theta and error for Euler method
    omega = 0.0
    theta = theta0
    err = [0.0] * nsteps
    for i in range(nsteps -1):
        t[i+1] = t[i]+dt
        # Euler method
        omegatmp = omega
        omega = omega -(g/l)*theta*dt
        theta = theta + omegatmp * dt
        # Error
        exact = theta0*math.cos(omega0*t[i+1])
        err[i+1] = abs(theta-exact)
        # Plotting the errors at this resolution
        pylab.loglog(t[refine::refine],err[refine::refine],'.-',label='dt = ' +str(dt))
pylab.legend(loc=4)
pylab.xlabel('Time - log scale')
pylab.ylabel('Absolute Error - log scale')
pylab.title('Euler method integration error')
pylab.grid(linestyle='-', which='major')
pylab.grid(which='minor')
pylab.show()







