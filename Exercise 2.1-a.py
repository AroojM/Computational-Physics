
# Exercise 2.1-a

import pylab

Nuclei_0 = input('Initial number of type A nuclei ->')
Tau_A = input('Decay time constant for type A->')
Tau_B = input('Decay time constant for type B->')
dt = input('Time step ->')
tmax = input('Time to end of simulation ->')

nsteps = int(tmax / dt)
Nuclei_A = [0.0] * nsteps
Nuclei_B = [0.0] * nsteps
t = [0.0] * nsteps

# Use Euler's method to integrate equations for radioactive decay

t[0] = 0.0 
Nuclei_A[0] = Nuclei_0
Nuclei_B[0] = 0.0

for i in range(nsteps-1):
    t[i+1] = t[i] + dt
    Nuclei_A[i+1] = Nuclei_A[i] - (Nuclei_A[i] / Tau_A) * dt
    Nuclei_B[i+1] = Nuclei_B[i] + ((Nuclei_A[i] / Tau_A)-(Nuclei_B[i] / Tau_B)) * dt

pylab.plot(t, Nuclei_A, 'o-r', label='Nuclei A - Tau_A = 1')
pylab.plot(t, Nuclei_B, 'o-b', label='Nuclei B - Tau_B = 10')
pylab.legend()
pylab.xlabel('Time')
pylab.ylabel('Nuclei A and Nuclei B')
pylab.title('Radioactive Decay')
pylab.grid()
pylab.show()
