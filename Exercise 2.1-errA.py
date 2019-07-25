
# Exercise 2.1 - Error_A

import pylab
import math
Nuclei_0 = input('Initial number of type A nuclei ->')
Tau_A = input('Decay time constant for type A->')
Tau_B = input('Decay time constant for type B->')
dtlow = input('Lowest resolution time step ->')
nres = input('Number of resolution refinements ->')
tmax = input('Time to end of simulation ->')
for n in range(nres):
    refine = 10**n
    dt = dtlow / refine
    nsteps = int(tmax / dt)
    Nuclei_A = Nuclei_0
    Nuclei_B = 0.0
    err_A = [0.0]*nsteps
    err_B = [0.0]*nsteps
    t = [0.0]*nsteps 

# Use Euler's method to integrate equations for radioactive decay
    for i in range(nsteps-1):
        t[i+1] = t[i] + dt
        Nuclei_A = Nuclei_A - (Nuclei_A / Tau_A) * dt
        exact_A = Nuclei_0 * math.exp(-t[i+1] / Tau_A)
        err_A[i+1] = abs((Nuclei_A - exact_A)/exact_A)
    #plot the error_A at this resolution
    pylab.loglog(t[refine::refine],err_A[refine::refine],
                 '.-',label='dt = ' + str(dt))
pylab.legend(loc=0)
pylab.xlabel('Time')
pylab.ylabel('Fractional error in A')
pylab.title('Radioactive decay Error_A')
pylab.grid(linestyle='-', which='major')
pylab.grid(which='minor')
pylab.show()


    
