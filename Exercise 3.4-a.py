# Exercise 3.4(a) - relax method

import math, pylab, mpl_toolkits.mplot3d

eps = 1e-2 # fractional error allowed
L = 1.0 # length of each side
N =  21 # Number of grid points on a side

dy = dx = L/(N-1.0)

xmax = ymax = 1.0
xmin = ymin = 0.0
# setting up arrays
x = pylab.array(range(N))*dx
y = pylab.array(range(N))*dy
x, y = pylab.meshgrid(x,y)
u0 = pylab.zeros((N, N))
u1 = pylab.zeros((N, N))
Ex = pylab.zeros((N,N))
Ey = pylab.zeros((N,N))

# boundary conditions
for j in range(N):
    u1[j, N-1] = u0[j, N-1] = u1[j, 0] = u0[j, 0] = 0.0
for k in range(N):
    u1[N-1, k] = u0[N-1, k] = u1[0, k] = u0[0, k] = 0.0
for j in range(5,16):
    u0[j,5] = -1.0
for k in range(5,16):
    u0[k,15] = 1.0

n = 0 # number of iterations
err = 1.0 # average error per site

while err > eps:
    # next iteration in refinement
    n = n+1
    err = 0.0
    for j in range(1, N-1):
        for k in range(1, N-1):
            u1[j,k] = (u0[j-1,k]+u0[j+1,k]+u0[j,k-1]+u0[j,k+1])/4.0
            err += abs(u1[j,k] - u0[j,k])
    err /= N**2
    u0, u1 = u1, u0 # swap old and new arrays for next iteration
    #internal boundary conditions
    for l in range(5,16):
        u0[l,5] = -1.0
    for m in range(5,16):
        u0[m,15] = 1.0

# calculate electric field
for i in range(1,N-1):
    for j in range(1,N-1):
        Ex[i,j] = -(u0[i+1,j]-u0[i-1,j])/ (2*dx)
        Ey[i,j] = -(u0[i,j+1]-u0[i,j-1])/ (2*dy)

# plot electric field
pylab.figure()
EF = pylab.quiver(y,x,Ex,Ey)
pylab.axis([xmin-dx, xmax+dx, ymin-dy, ymax+dy])
pylab.title('Electric Field - relax method')
pylab.xlabel('x')
pylab.ylabel('y')
pylab.show()
        
## surface plot of final solution
fig = pylab.figure()
axis = fig.gca(projection = '3d', azim= -60, elev = 20)
surf = axis.plot_surface(x,y, u0.T, rstride=1, cstride=1, cmap=pylab.cm.jet)
axis.contour(x, y, u0.T, 10, zdir='z', offset=-1.0)
axis.set_xlabel('x')
axis.set_ylabel('y')
axis.set_zlabel('u')
axis.set_zlim(-1.0,1.0)
fig.colorbar(surf)
pylab.title('Electric Potential - relax method')
pylab.show()
