# Exercise 4.1

import pylab, math, random

L1 = 100 #length of horizontal side of column
L2 = 10000 #length of vertical side of column

M1 = 40 #horizontal length of rectangle of particles
M2 = 20 #vertical length of rectangle of particles

B1 = 10 #number of coarse-grating bins per horizontal side
B2 = 10 #number of coarse-grating bins per vertical side


nsteps = input('Number of steps in walk ->')
steps = range(nsteps)
seed = input('Random number seed ->')
random.seed(seed)

#intial positions of particles form a M1*M2 block in the middle of the column
xside = range((L1-M1)//2,(L1+M1)//2)
yside = range((L2-M2)//2,(L2+M2)//2)
x = [i for i in xside for j in yside] # x-locations of the particles
y = [j for i in xside for j in yside] # y-locations of the particles
N = len(xside) * len(yside) # number of particles
S = [0.0]*nsteps # entropy
Num1 = pylab.zeros((B1,B2)) # number of particle  in each bin

#setup animated figure
pylab.figure(figsize=(5,10))
points, = pylab.plot(x,y,',')
pylab.xlim(0,L1)
pylab.ylim(0,L2)
pylab.xticks(range(0,L1+B1,B1))
pylab.yticks(range(0,L2+B2,1000))
pylab.xlabel('x - (m)')
pylab.ylabel('y - (m)')
pylab.grid()
Z = [0.0]*nsteps # partition function
H = [0.0]*B2 # list for height (at the middle of each bin B2)

# simulate the random walks
for n in steps:
    # update plots
    points.set_data(x,y)
    pylab.title('step %d' % n)
    pylab.pause(0.001)
    #pylab.draw()

    # update positions of particles and update counts in bins
    Num1.fill(0)
    for i in range(N):
        dx = random.choice([-1,1])
        # choosing random number for biased walk in vertical direction
        if random.random() < 0.475:
            dy = 100
        else:
            dy = -100
        x[i] += dx
        y[i] += dy
        # make sure that the particles stay in the column
        if x[i] < 0 or x[i] >= 100: x[i] -= dx
        if y[i] < 0 or y[i] >= 10000: y[i] -= dy
        # increment count in bin containing particle
        Num1[x[i]//10, y[i]//1000] += 1.0
    # number of particles with same energy (i.e. in bin with same B2 value)
    Num2 = [0.0]*B2
    for s in range(B2):
        for r in range(B1):
            Num2[s] += Num1[r,s]
    for s in range(B2):
        H[s] = (s+1)*1000 - 500
        Z[n] += math.exp(- 0.000125*H[s])*Num2[s]
for n in steps:
    S[n] = math.log(Z[n])

pylab.figure()
pylab.plot(steps, S)
pylab.title('Entropy versus Steps')
pylab.xlabel('Step')
pylab.ylabel('Entropy')
pylab.show()

pylab.figure()
pylab.plot(H, Num2, 'o')
pylab.title('Number of partilces versus height')
pylab.xlabel('Height - (m)')
pylab.ylabel('Number of particles')
pylab.show()

 

