# Ex 3.2 - implicit scheme - graphs

import math
import pylab
# parameters
D = 1.0 #thermal diffusivity, m^2/s
L = 1.0 # length of rod, m
dt = 0.0001 #time step
tmax = 0.02 #evolution time
N = 100 #number of grid points
dx = L/(N) # grid size
#coefficients of the tridiagonal matrix
alpha = gamma = -D * dt / dx**2
beta = 1.0 + 2.0 * D * dt / dx**2
#creating lists
x = [0.0]*(N+1)
v = [0.0]*(N+1)
u = [0.0]*(N+1)
#intial condition - delta function is in center
u[50] = 1/dx
#spatial location of grid points
for j in range(N+1):
    x[j] = (j * dx) - (L / 2)
#perform the evolution, implicit differencing scheme
t = 0.0
p = 0
tstep = 0.0040
pylab.figure()
while t < tmax:
    #swap u and v
    u, v = v, u
    #set the j=1 and j=N-1 points of v to the correct values
    v[1] -= alpha * u[0]
    v[N-1] -= gamma * u[N]
    #forward sweep
    u[1] = v[1] / beta
    v[1] = gamma / beta
    for j in range(2,N):
        den = beta - alpha*v[j-1]
        u[j] = (v[j] - alpha*u[j-1])/den
        v[j] = gamma / den
    #backward sweep
    for j in reversed(range(1,N-1)):
        u[j] -= u[j+1] * v[j]
    t += dt
    if t >= p*tstep:
        pylab.plot(x, [z  for z in u], label='t (s)='+str(round(t,4)) )
        pylab.xlim(-L/2,L/2)
        pylab.xlabel('x(m)')
        pylab.ylabel('Temperature (Celcius)')
        pylab.legend(loc=(1.03,0.2))
        p += 1

pylab.show()
        












    

