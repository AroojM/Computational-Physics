# Exercise 2.2-b (Launch angle for greatest range
#                 for projectile WITH air drag)

import math
import pylab
g = 9.8 # acceleration due to gravity (m/s^2)
v0 = 10.0 # initial velocity (m/s)
b = 0.04 # drag constant (m^-1)
angles = [0.0]*91 # create list for launch angles (degrees)
angles[0] = 0.0
for k in range(91):
    angles[k] = angles[0] + 1.0*k
dt = 0.01 # time step (s)
xmax = 0.0 # greatest range
thetamax = 0.0 # angle for greatest range
for theta in angles:    
    x = [0.0]
    y = [0.0]
    vx = [v0 * math.cos(theta*math.pi/180.0)]
    vy = [v0 * math.sin(theta*math.pi/180.0)]
    v = [v0]
    # using Euler's method to integrate projectile equations of motion
    i = 0
    while y[i]>= 0.0 :
        #extend the lists
        x.append(0.0)
        y.append(0.0)
        vx.append(0.0)
        vy.append(0.0)
        v.append(0.0)
        x[i+1] = x[i] + vx[i]*dt
        y[i+1] =y[i] + vy[i]*dt
        vx[i+1] = vx[i] - b*v[i]*vx[i]*dt
        vy[i+1] = vy[i] - g*dt - b*v[i]*vy[i]*dt
        v[i+1] = math.sqrt((vx[i+1])**2+(vy[i+1])**2)
        i = i+1
    # find location of point where projectile hits ground  
    xground = x[i-1] - y[i-1]*(x[i]-x[i-1])/(y[i]-y[i-1])
    if xground > xmax:
        xmax = xground
        thetamax = theta
print "For initial speed", v0, "m/s and drag constant", b, "m^-1"
print "Greatest range achieved is = ", round(xmax,2), "m"
print "The launch angle for this greatest range is =", thetamax, "degrees"



        
