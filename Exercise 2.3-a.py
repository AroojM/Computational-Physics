# Exercise 2.3-a

import pylab
g = 9.8 # acceleration due to gravity (m/s^2)
l = input("Length of pendulum (meters) ->")
m = input("Mass of the pendulum bob (kg) ->")
theta0 = input("Initial angle (radians) ->")
dt = input("Time step (seconds) ->")
tmax = input("Time to end of simulation (seconds) ->")
nsteps = int(tmax/dt)
t = [0.0] * nsteps
# omega, theta and energy for Euler method
omega1 = [0.0] * nsteps
theta1 = [0.0] * nsteps
energy1 = [0.0] * nsteps
# omega, theta and energy for Euler-Cromer method
omega2 = [0.0] * nsteps
theta2 = [0.0] * nsteps
energy2 = [0.0] * nsteps

theta1[0] = theta0
theta2[0] = theta0
for i in range(nsteps - 1):
    # Euler method 
    omega1[i+1] = omega1[i] -(g/l)*theta1[i]*dt
    theta1[i+1] = theta1[i] + omega1[i] * dt
    #  Euler-Cromer method 
    omega2[i+1] = omega2[i] -(g/l)*theta2[i]*dt
    theta2[i+1] = theta2[i] + omega2[i+1] * dt
    t[i+1] = t[i] + dt
# Compute energy at each step for Euler and Euler-Cromer method
for k in range(nsteps):
    energy1[k] =(m*l/2.0)*(l*omega1[k]**2 + g*theta1[k]**2)
    energy2[k] =(m*l/2.0)*(l*omega2[k]**2 + g*theta2[k]**2)
# Plotting energy vs time
pylab.figure()
pylab.title('Total energy simple pendulum - Euler-Cromer method')
pylab.plot(t, energy1, label='Euler')
pylab.plot(t, energy2, label='Euler-Cromer')
pylab.xlabel('Time (s)')
pylab.ylabel('Energy (J)')
pylab.legend(loc=9)
pylab.show()







