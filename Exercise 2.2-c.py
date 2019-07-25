
# Exercise 2.2-c (Launch angle for greatest range VS Drag constant)

import math
import pylab
g = 9.8 # acceleration due to gravity (m/s^2)
v0 = 10.0 # initial velocity (m/s)
drag=[0.0, 0.01, 0.02, 0.03,0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1] # drag constant (m^-1)
angles = [0.0]*182 # create list for launch angles (degrees)
angles[0] = 0.0
for k in range(182):
    angles[k] = angles[0] + 0.5*k
dt = 0.01 # time step (s)
thetamax = [0.0]*11 # create list for thetamax (angle for greatest range)
j = 0
print "b (m^-1)","Theta max(deg)","  Range(m)"
for b in drag:
    xmax = 0.0 # greatest range
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
            thetamax[j] = theta        
    print b,"        ",thetamax[j], "       ",round(xmax,2)
    j = j + 1
pylab.scatter(drag, thetamax, s=8, facecolors='r', edgecolors='r')
pylab.title('Launch angle for maximum range VS Drag Constant')
pylab.xlabel('Drag Constant - b (m^-1)')
pylab.ylabel('Thetamax (degrees)')
pylab.axis([-0.01,0.11,36.0,49.0])
pylab.show()




