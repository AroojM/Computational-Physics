# Exercise 2.7

import math, pylab

def rk45(t, x, f, dt):
    """
    Fifth order and embedded fourth order Runge-Kutta integration step with
    error estimate
    """
    # Cash - Karp parameters
    a2, a3, a4, a5, a6 = 1./5., 3./10., 3./5., 1., 7./8.
    b21 = 1./5.
    b31, b32 = 3./40., 9./40.
    b41, b42, b43 = 3./10., -9./10., 6./5.
    b51, b52, b53, b54 = -11./54., 5./2., -70./27., 35./27.
    b61, b62, b63, b64, b65 = 1631./55296., 175./512., 575./13824., 44275./110592., 253./4096.
    c1, c2, c3, c4, c5, c6 = 37./378., 0., 250./621., 125./594., 0., 512./1771.
    d1, d2, d3, d4, d5, d6 = 2825./27648., 0., 18575./48384., 13525./55296., 277./14336., 1./4.
    e1, e2, e3, e4, e5, e6 = c1-d1, c2-d2, c3-d3, c4-d4, c5-d5, c6-d6

    # evaluate the function at the six points
    dx1 = f(t,x) * dt
    dx2 = f(t + a2*dt, x + b21*dx1) * dt
    dx3 = f(t + a3*dt, x + b31*dx1 + b32*dx2) * dt
    dx4 = f(t + a4*dt, x + b41*dx1 + b42*dx2 + b43*dx3 ) * dt
    dx5 = f(t + a5*dt, x + b51*dx1 + b52*dx2 + b53*dx3 + b54*dx4) * dt
    dx6 = f(t + a6*dt, x + b61*dx1 + b62*dx2 + b63*dx3 + b64*dx4 + b65*dx5) * dt
    # compute and return the error and the new value of x
    err = e1*dx1 + e2*dx2 + e3*dx3 + e4*dx4 + e5*dx5 + e6*dx6
    return x + c1*dx1 + c2*dx2 + c3*dx3 + c4*dx4 + c5*dx5 + c6*dx6, err

def ark45(t, x, f, dt, epsabs = 1e-6, epsrel = 1e-6):
    """
    Adaptive Runge-Kutta integration step.
    """

    safe = 0.9 # safety factor for step estimate
    # compute the required error
    e0 = epsabs + epsrel * max(abs(x))
    dtnext = dt
    while True:
        # take a step and estimate the error
        dt = dtnext
        result, error = rk45(t, x, f, dt)
        e = max(abs(error))
        dtnext = dt * safe * (e0 / e)**0.2
        if e < e0: # accept step: return x, t and dt for next step
            return result, t + dt, dtnext

# defining a function that returns F for equations of motion of particles
def F(t,X):
    x1 = X[0]
    vx1 = X[1]
    y1 = X[2]
    vy1 = X[3]
    x2 = X[4]
    vx2 = X[5]
    y2 = X[6]
    vy2 = X[7]
    x3 = X[8]
    vx3 = X[9]
    y3 = X[10]
    vy3 = X[11]
    r12 = math.sqrt((x1-x2)**2+(y1-y2)**2) # distance between particle 1 & 2
    r13 = math.sqrt((x1-x3)**2+(y1-y3)**2) # distance between particle 1 & 3
    r23 = math.sqrt((x2-x3)**2+(y2-y3)**2) # distance between particle 1 & 4
    ax1 = -(m2 *(x1-x2) / r12**3)-(m3 * (x1-x3) / r13**3)
    ay1 = -(m2 *(y1-y2) / r12**3)-(m3 * (y1-y3) / r13**3)
    ax2 = -(m1 *(x2-x1) / r12**3)-(m3 * (x2-x3) / r23**3)
    ay2 = -(m1 *(y2-y1) / r12**3)-(m3 * (y2-y3) / r23**3)
    ax3 = -(m2 *(x3-x2) / r23**3)-(m1 * (x3-x1) / r13**3)
    ay3 = -(m2 *(y3-y2) / r23**3)-(m1 * (y3-y1) / r13**3)
    return pylab.array([vx1, ax1, vy1, ay1, vx2, ax2, vy2, ay2, vx3, ax3, vy3, ay3])

#data
# masses
m1 = 3
m2 = 4
m3 = 5

# initial conditions
# particle 1
x10 = 4.0
y10 = 0.0
vx10 = 0.0
vy10 = 0.0
# particle 2
x20 = 0.0
y20 = 3.0
vx20 = 0.0
vy20 = 0.0
# particle 3
x30 = 0.0
y30 = 0.0
vx30 = 0.0
vy30 = 0.0

x1 = [x10]
y1 = [y10]
x2 = [x20]
y2 = [y20]
x3 = [x30]
y3 = [y30]

dt =0.01
tmax = 100.0
t = [0.0]

# for plotting output at different intervals
p = 1
ptime = 0.0
pstep = 10

# integrate Newton's equatins of motion using ark45
# X is a vector that contains the positions and velocities being integrated
X = pylab.array([x10, vx10, y10, vy10, x20, vx20, y20, vy20, x30, vx30, y30, vy30])
T = 0.0
while T < tmax:
    X, T, dt = ark45(T, X, F, dt)
    x1.append(X[0])
    y1.append(X[2])
    x2.append(X[4])
    y2.append(X[6])
    x3.append(X[8])
    y3.append(X[10])
    t.append(T)
    if ptime >= p*pstep and ptime < tmax:
        pylab.figure()
        pylab.plot(x1, y1, 'k')
        pylab.plot(x2, y2, 'b')
        pylab.plot(x3, y3, 'r')
        pylab.title('Trajectories of three particles at time = %.2f ' % ptime)
        pylab.xlabel('x')
        pylab.ylabel('y')
        pylab.grid()
        p +=1
    ptime += dt

pylab.show()






    
