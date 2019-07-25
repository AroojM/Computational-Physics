# Exercise 2.6 - b

import math, pylab
GM = (2.0 * math.pi)**2 #Heliocentric gravitational constant
A = [0.0002,0.0003, 0.0005, 0.001, 0.002,0.003] # Constant in general relativity force law
# Defining a function that implements RK4 integration
def RK4(t,z,g,dt):
    k1 = g(t,z) * dt
    k2 = g(t+0.5*dt, z+0.5*k1) * dt
    k3 = g(t+0.5*dt, z+0.5*k2) * dt
    k4 = g(t+dt, z+k3) * dt
    return z + k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0
# Defining a function that returns the vector F(t,X) in
# equations of motion for planetary orbits
def F(t,X):
    x = X[0]
    vx = X[1]
    y = X[2]
    vy = X[3]
    r = math.sqrt(x**2 + y**2)
    ax = -GM * x * ( 1.0/ r**3 + alpha / r**5)
    ay = -GM * y * ( 1.0/ r**3 + alpha / r**5) 
    return pylab.array([vx, ax, vy, ay])

x0 = 0.47 #Initial x position (AU)
y0 = 0.0 #Initial y position (AU)
vx0 = 0.0 #Initial x velocity (AU/yr)
vy0 = 8.2 #Initial y velocity (AU/yr)
dt = 0.0001 #Time step (yr)
tmax = 2.0 #Time to the end of simulation (yr)
rate = [0.0] * 6 # rate of perihelion advance - d(theta)/dt
#- list of 6 because 6 values of alpha
p = 0 # index for rate list
# integrating equations for different values of alpha
for alpha in A:
    nsteps = int(tmax / dt)
    x = [0.0] * nsteps
    y = [0.0] * nsteps
    R = [0.0] * nsteps # distance of mercury from sun
    sR = [0.0] * (nsteps-1) # rate of change of R
    theta = [0.0] # angle at which sR changes sign
    time = [0.0] # time at which sR changes sign
    # Integrating General Relativity equations of motion using RK4
    # X ia a vector that contains the positions and velocities being integrated
    X = pylab.array([x0, vx0, y0, vy0])
    for i in range(nsteps):
        x[i] = X[0]
        y[i] = X[2]
        R[i] =  math.sqrt(x[i]**2 + y[i]**2)
        # Updating the vector X to next time step    
        X = RK4(i*dt, X, F, dt)
     # slope of R
    for j in range(1,nsteps-1):
        sR[j] = (R[j+1] - R[j-1])/dt
    # determine R, x and theta at point where slope of R changes sign
    for k in range(nsteps -1):
        if sR[k]<0 and  sR[k-1] > 0:
            nexttheta = math.acos(x[k] / R[k])  
            theta_deg = nexttheta*180.0 /(math.pi)  # convert angle to degrees
            theta.append(theta_deg) # extend list
            nexttime = dt*k
            time.append(nexttime)  # extend list
    # start - linear fit for theta vs time
    c = time
    d = theta
    cSQ = [ i**2 for i in c]
    cd = [i*j for i,j in zip(c,d)]
    cbar = math.fsum(c) / len(c)
    dbar = math.fsum(d) / len(d)
    cSQbar = math.fsum(cSQ) / len(cSQ)
    cdbar = math.fsum(cd) / len(cd)
    cbarSQ = cbar **2
    M = (cdbar -cbar*dbar)/(cSQbar - cbarSQ)
    C = ((cSQbar * dbar)-(cbar*cdbar))/(cSQbar - cbarSQ)
    # end - linear fit for theta vs time
    
    rate[p] = M
    print 'Rate of perihelion advance for alpha = ',alpha, 'is', rate[p], 'degrees/yr'
    p = p + 1

# start - linear fit for alpha vs rate
c = A
d = rate
cSQ = [ i**2 for i in c]
cd = [i*j for i,j in zip(c,d)]
cbar = math.fsum(c) / len(c)
dbar = math.fsum(d) / len(d)
cSQbar = math.fsum(cSQ) / len(cSQ)
cdbar = math.fsum(cd) / len(cd)
cbarSQ = cbar **2
M = (cdbar -cbar*dbar)/(cSQbar - cbarSQ)
C = ((cSQbar * dbar)-(cbar*cdbar))/(cSQbar - cbarSQ)
# end - linear fit for alpha vs rate

# slope of d(theta)/dt VS alpha
print 'Slope of Rate of perihelion advance VS Alpha is', M , '(degrees per year per unit alpha)'
# extrapolate the slope to mercury value
MercuryRate1 = M * 1.1 * 1e-8
MercuryRate2 = MercuryRate1 * 100.0 * 3600.0
print 'Rate of perihelion advance for Mercury is', MercuryRate1, 'degrees/year', '=', round(MercuryRate2,3), 'arcseconds/century'

pylab.figure()
pylab.plot(A, rate,'o')
pylab.title('Rate of perihelion advance VS Alpha')
pylab.xlabel('alpha')
pylab.ylabel('rate of perihelion advance - (degrees/yr)')
pylab.grid()
pylab.show()
