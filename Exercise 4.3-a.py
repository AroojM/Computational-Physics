# Exercise 4.3 - a - V-shaped potential

import math, random, pylab
# parameters
hbar = 1.0
m = 1.0
# input number of grid points, steps, random seed
##N = input('Number of grid points ->')
##nstep = input('Number of steps ->')
##seed =  input('Random number seed ->')
N = 21
nstep = 500
seed = 35
random.seed(seed)
Step = [0.0]*(nstep+1) # list for plotting E
j = 0
while j < nstep:
    Step[j+1] = j+1
    j += 1
# setup grid and initial guess
xmin = -5.0
xmax = +5.0
dx = (xmax - xmin)/(N-1)
x = pylab.arange(xmin, xmax + 0.1*dx, dx)
psi = pylab.ones(N) #initial guess
psi[0] = psi[N-1] = 0.0 #end points fixed at zero
# compute energy, potential, normalization
V = pylab.zeros(N)
E = pylab.zeros(nstep + 1)
ssq = 0
for i in range(1,N-1):
    V[i] = abs(x[i]) # V-shaped potential
    H = -hbar**2*(psi[i-1]-2.0*psi[i]+psi[i+1])/(2*m*dx**2)
    H += V[i]*psi[i]
    E[0] += psi[i]*H*dx
    ssq += psi[i]**2 *dx
E[0] /= ssq
psi /= ssq**0.5
# plot
pylab.figure()
n=1
p = 0
pstep = 100
while n<= nstep:
    if n >= p*pstep:
        pylab.plot(x, [u  for u in psi],'o-', label='step ='+str(p*pstep))
        pylab.title('Wave function - V-shaped potential')
        pylab.xlabel('x')
        pylab.ylabel('$\psi$')
        pylab.legend(loc=(1.03,0.2))
        p += 1
    # choose a random point and a random amount to change psi
    tmp = pylab.copy(psi) #temporary wavefunction trial
    j = random.choice(range(1,N-1))
    tmp[j] *= random.uniform(0.8, 1.2)
    # normalize and compute energy
    E[n] = 0.0
    ssq = 0.0
    for i in range(1, N-1):
        H = -hbar**2*(tmp[i-1] - 2.0*tmp[i] + tmp[i+1])/(2*m*dx**2)
        H += V[i]*tmp[i]
        E[n] += tmp[i]*H*dx
        ssq += tmp[i]**2 *dx
    E[n] /= ssq
    # test if the trial wavefunction reduces energy
    if E[n] < E[n-1]:
        # update current wavefunction
        psi = tmp / ssq**0.5 
        #increment step count
        n += 1
pylab.figure()
pylab.plot(Step,E)
pylab.title('Energy - V-shaped potential')
pylab.ylabel('Energy ')
pylab.xlabel('Step number')
pylab.grid()
pylab.show()
    
