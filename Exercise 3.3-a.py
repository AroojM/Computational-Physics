# Exercise 3.3-a

import math, cmath, pylab
hbar = 1.0 # reduced Plank's constant
m = 1.0 # mass
k0 = 1.0 # initial wavenumber
p0 = hbar * k0 # initial momentum
A = 50.0
# grid and time intervals
dx = 0.25/k0
dt = 0.25*m/(hbar*k0**2)
tmax = 150.0*m/(hbar*k0**2)
xmin = -200.0/k0
xmax = 250.0/k0
N = int((xmax - xmin)/dx)
# initial data
x = [0.0]*(N+1)
u = [0.0]*(N+1)
v = [0.0]*(N+1)
p = [0.0]*(N+1)
x0 = -100.0/k0
sigma = 10.0/k0
for j in range(N+1):
    x[j] = xmin + j*dx
    u[j] = cmath.exp(-(x[j]-x0)**2/(4.0*sigma**2)+1j*k0*x[j])
    u[j] = u[j]/(2.0*math.pi*sigma**2)**0.25
    p[j] = u[j].real**2 + u[j].imag**2
# potential
E0 = -A*(hbar*k0)**2 / (2.0*m)
V = [0.0]*(N+1)
for j in range(N+1):
    if x[j] > 0: V[j] = E0
# setup coefficients of the tridiagonal matrix
alpha = gamma = -1j*hbar*dt/(4.0*m*dx**2)
beta = [0.0]*(N+1)
for j in range(N):
    beta[j] = 1.0 -2.0*alpha+1j*(V[j]/(2.0*hbar))*dt
# prepare plots
fig = pylab.figure()
pylab.title('Scattering of wave packet from potential shelf in time increments of 50 seconds, Alpha = %0.1f' % A)
ax1 =fig.add_subplot(111)
ax1.set_xlim(xmin, xmax)
ax1.set_ylim(-0.035, 4.0*max(p))
ax1.set_xlabel('x')
ax1.set_ylabel('Probability density')
ax2 = ax1.twinx()
ax2.set_xlim(xmin, xmax)
ax2.set_ylim(min(V)-0.05, -7.0*min(V))
ax2.set_ylabel('V')
ax2.plot(x, [Vj for Vj in V],'k')
# perform the evolution
t = 0.0
n = 0
tstep = 50
while t  < tmax + 0.5*dt:
    for j in range(N+1):
        p[j] = u[j].real**2 + u[j].imag**2
    #plotting 
    if t >= n*tstep:
        ax1.plot(x, [P+n*0.042 for P in p])
        n += 1
    # set the values of the rhs
    for j in range(1,N):
        v[j] = -alpha*u[j-1]+(2.0-beta[j])*u[j]-gamma*u[j+1]
    v[1] -= alpha*u[0]
    v[N-1] -= gamma*u[N]
     # forward sweep
    u[1] = v[1]/beta[1]
    v[1] = gamma/beta[1]
    for j in range(2,N):
        den =  beta[j] - alpha*v[j-1]
        u[j] = (v[j] - alpha*u[j-1])/den
        v[j] = gamma/den
    # backward sweep
    for j in reversed(range(1,N)):
        u[j] -= u[j+1]*v[j]
    t += dt
pylab.show()
















    
