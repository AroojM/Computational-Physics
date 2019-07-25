# Exercise 2.4 - b - chaotic motion

import math
import pylab

g = 9.8 #acceleration due to gravity (m/s^2)
l = 1.0
m = 1.0
theta1= [1.57,1.58, 2.8, 3.8, 4.7, 4.8]
theta20 = 0.0
dt = 0.0001
tmax = 20.0

nsteps =  int(tmax / dt)

for theta10 in theta1:
    t = [0.0] * nsteps
    p1 = [0.0] * nsteps
    p2 = [0.0] * nsteps
    q1 = [0.0] * nsteps
    q2 = [0.0] * nsteps
    x1 = [0.0] * nsteps
    x2 = [0.0] * nsteps
    y1 = [0.0] * nsteps
    y2 = [0.0] * nsteps

#initialize
    q1[0] = theta10
    q2[0] = theta20
    x1[0] = l * math.sin(q1[0])
    y1[0] = -l * math.cos(q1[0])
    x2[0] = l * (math.sin(q1[0])+math.sin(q2[0]))
    y2[0] = -l * (math.cos(q1[0])+math.cos(q2[0]))

# Using Euler-Cromer method to integrate the equations for double pendulum
    for i in range(nsteps - 1):
        s = math.sin(q1[i] - q2[i])
        c = math.cos(q1[i] - q2[i])
        D = m*(l**2)*(1 + s**2)
        A = p1[i] * p2[i] * s / D
        B = m*(l**2)*s*c*(p1[i]**2 + 2*p2[i]**2 - p1[i]*p2[i]*c) / (D**2)
        p1[i+1] = p1[i]+ (-2*m*g*l*math.sin(q1[i]) - A + B)*dt
        p2[i+1] = p2[i] +(-m*g*l*math.sin(q2[i]) + A - B)*dt
        q1[i+1] = q1[i] + ((p1[i+1]-p2[i+1]*c)/ D)*dt
        q2[i+1] = q2[i] + ((2*p2[i+1]-p1[i+1]*c)/ D)*dt
        t[i+1] =t[i] + dt

        #put q1 and q2 in range -pi to +pi
        q1[i+1] = (q1[i+1]+math.pi) % (2.0 * math.pi) - math.pi
        q2[i+1] = (q2[i+1]+math.pi) % (2.0 * math.pi) - math.pi

        #Computing (x,y): locations of pendulum 1 and 2
        x1[i+1] = l * math.sin(q1[i+1])
        y1[i+1] = -l * math.cos(q1[i+1])
        x2[i+1] = l * (math.sin(q1[i+1])+math.sin(q2[i+1]))
        y2[i+1] = -l * (math.cos(q1[i+1])+math.cos(q2[i+1]))

    label = theta10 
     #Plotting results
    pylab.figure()
    pylab.title('Double pendulum -theta10 = %f'% theta10)
    pylab.plot(t, q2, label='theta2')
    pylab.plot(t, q1, label='theta1')
    pylab.xlabel('t(s)')
    pylab.ylabel('Angle(rad)')
    pylab.legend(loc=9)
    pylab.grid()

    pylab.figure(figsize=(6,6))
    pylab.title('Lissajous curves - theta10 = %f'% theta10) 
    pylab.plot(q1, q2)
    pylab.xlabel('theta1 (rad)')
    pylab.ylabel('theta2 (rad)')
    minmax = max(abs(min(q1+q2)),abs(max(q1+q2)))
    pylab.axis([-minmax,minmax,-minmax,minmax],aspect='equal')
    pylab.grid()
    
    pylab.figure()
    pylab.title('Double pendulum trace - theta10 = %f ' % theta10)
    pylab.plot(x2, y2)
                
    pylab.xlabel('x (m)')
    pylab.ylabel('y (m)')
##    pylab.axis([-2.1,2.1,-2.1,2.1],aspect='equal')
    pylab.grid()

pylab.show()


























