# Exercise 2.2-a (Trajectory of projectile WITH air drag)

import math
import pylab
g = 9.8 # acceleration due to gravity (m/s^2)
v0 = 10.0 # initial velocity (m/s)
b = 0.04 # drag constant (m^-1)
angles = [30.0, 35.0, 40.0, 45.0, 50.0, 55.0] # launch angles (degrees)
dt = 0.01 # time step (s)
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
    pylab.plot(x, y, label=str(theta)+' degrees')
    pylab.title('Trajectory of a projectile with air drag')
pylab.xlabel('x(m)')
pylab.ylabel('y(m)')
pylab.ylim(ymin=0.0)
pylab.axis([0.0,11.0,0.0,4.5])
pylab.legend()
pylab.show()
        
