# Exercise 2.6 - a

import math, pylab
GM = (2.0 * math.pi)**2 #Heliocentric gravitational constant
alpha = 0.01 # Constant in general relativity force law
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
tmax = 1.5 #Time to the end of simulation (yr)

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
        theta_deg = nexttheta*180.0 /(math.pi) # convert angle to degrees
        theta.append(theta_deg) # extend list
        nexttime = dt*k
        time.append(nexttime) # extend list
##print theta
pylab.figure()
pylab.plot(time, theta, 'o')
pylab.title('T vs theta for alpha %3f' %alpha)
pylab.xlabel('time - yr')
pylab.ylabel('theta - degrees')

pylab.figure(figsize=(6,6))
pylab.plot(x, y, '-')
pylab.title('Precession of perihelion - alpha = 0.01')
pylab.xlabel('x (AU)')
pylab.ylabel('y (AU)')
minmax = 1.1 * max(abs(min(x+y)),abs(max(x+y)))
pylab.axis([-minmax,minmax,-minmax,minmax], aspect='equal')
pylab.grid()
pylab.show()
