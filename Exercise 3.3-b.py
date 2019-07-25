# Exercise 3.3-b

import math, cmath, pylab
hbar = 1.0 # reduced Plank's constant
m = 1.0 # mass
k0 = 1.0 # initial wavenumber
p0 = hbar * k0 # initial momentum
A = 70.0
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
# perform the evolution
t = 0.0
n = 0
tstep = 50
while t  < tmax + 0.5*dt:
    for j in range(N+1):
        p[j] = u[j].real**2 + u[j].imag**2
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
Prob1 = 0.0 # Probability that particle is reflected
Prob2 = 0.0 # Probability that particle is transmitted
for j in range(N+1):
    if x[j] < 0.0:
        Prob1 += p[j]*dx
    if 0.0 < x[j]:
        Prob2 +=p[j]*dx
print 'Probability that particle is reflected for Alpha =', A, 'is', round(Prob1,4)
##print Prob2
##print Prob1+Prob2
##print Prob1/(Prob1+Prob2)


















    
