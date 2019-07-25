# Ex 3.2 - ftcs
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
N = 60
dx = L/(N) # grid size
#creating lists
x = [0.0]*(N+1)
u0 = [0.0]*(N+1)
u1 = [0.0]*(N+1)
#intial condition - delta function is in center
u0[30] = 1/dx
#spatial location of grid points
for j in range(N+1):
    x[j] = (j * dx) - (L / 2)
#prepare animated plot
pylab.ion
line, = pylab.plot(x, u0, '-k')
pylab.xlim(-L/2,L/2)
pylab.xlabel('x(m)')
pylab.ylabel('Temperature (Celcius)')
#perform the evolution
t = 0.0
while t < tmax:
    #update plot
    line.set_ydata(u0)
    pylab.title('t = %5f' % t)
    pylab.draw()
    pylab.pause(0.1)
    #derivatives at interior points
    for j in range(1,N):
        u1[j] = u1[j] = u0[j] + dt*D*(u0[j+1] - 2.0*u0[j] +u0[j-1])/dx**2
    #swap old and new lists
    u0, u1 = u1, u0
    t += dt
#freeze final plot
pylab.ioff()
pylab.show()
#   u1[j] = (1/(math.sqrt(4*math.pi*D*t)))*(math.exp(-(x[j])**2/(4*D*t)))
        












    

