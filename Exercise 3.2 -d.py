# Ex 3.2 - implicit
import math
import pylab
#fixed parameters
D = 1.0 #thermal diffusivity, m^2/s
#input parameters
#dx = input('Grid spacing->')
##dt = input('Time step->')
##tmax = input('Evolution time ->')
##N = input('Number of grid points->') #grid points
##dx = L/N # grid size
dt = 0.0001
tmax = 0.005
L = 1.0
N = 100
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
#prepare animated plot
pylab.ion
line, = pylab.plot(x, u, '-k')
pylab.xlim(-L/2,L/2)
pylab.xlabel('x(m)')
pylab.ylabel('Temperature (Celcius)')
#perform the evolution
t = 0.0
while t < tmax:
    #update plot
    line.set_ydata(u)
    pylab.title('t = %5f' % t)
    pylab.draw()
    pylab.pause(0.1)
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
#freeze final plot
pylab.ioff()
pylab.show()
#   u1[j] = (1/(math.sqrt(4*math.pi*D*t)))*(math.exp(-(x[j])**2/(4*D*t)))
        












    

