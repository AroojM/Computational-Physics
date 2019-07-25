# Exercise 4.2

import math, pylab, random
J = 1.0 # exchange energy
T = input('Temperature in unit of kT ->')
L = input('Number of atoms per side of lattice ->')
nsweep = input('Number sweeps to average ->')
seed = input('Random number seed ->')
random.seed(seed)
N = L**2
H = pylab.arange(-5.0, 5.0, 0.2)
e = pylab.zeros(len(H))
m = pylab.zeros(len(H))
#initial data
s = pylab.ones((L,L))
E = 0.0
M = 0.0

for h in range(len(H)):
    for i in range(L):
        for j in range(L):
            E +=-J*s[i,j]*(s[(i+1)%L,j]+s[i,(j+1)%L])- H[h]*s[i,j]
for i in range(L):
        for j in range(L):
            M += s[i,j]
for h in range(len(H)):
    #average nsweep sweeps
    for sweep in range(nsweep):
        #sweep over all particles in lattice
        for i in range(L):
            for j in range(L):
                #compute energy required to flip spin
                dE = 2.0*J*s[i,j]*(s[(i+1)%L,j]+s[(i-1)%L,j]+s[i,(j+1)%L]+s[i,(j-1)%L]) + 2*H[h]*s[i,j]
                #Metropolis algorithm to see if we should accept trial
                if dE <= 0.0 or random.random() <= math.exp(-dE/T):
                    #accept trial: reverse spin; return dE and dM
                    s[i,j] *= -1
                    M += 2.0 *s[i,j]
                    E += dE
        #update running means 
        deltam = M -m[h]
        m[h] += deltam / (sweep + 1)               
    m[h] /= N
#produce plots
pylab.figure()
pylab.plot(H, m, 'o-')
pylab.title('T = %d ' % T)
pylab.axis([-5.0,5.0,-1.2,1.2])
pylab.xlabel('H')
pylab.ylabel('M')
pylab.grid()
pylab.show()
















            
