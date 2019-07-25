
# Exercise 2.1-b

import pylab
import math
Nuclei_0 = input('Initial number of type A nuclei ->')
Tau_A = input('Decay time constant for type A->')
Tau_B = input('Decay time constant for type B->')
dt = input('Time step ->')
tmax = input('Time to end of simulation ->')
nsteps = int(tmax / dt)
Nuclei_A = Nuclei_0
Nuclei_B = 0.0
t = [0.0] * nsteps
err_A = [0.0]*nsteps
err_B = [0.0]*nsteps
# Use Euler's method to integrate equations for radioactive decay
t[0] = 0.0 
for i in range(nsteps-1):
    t[i+1] = t[i] + dt
    Nuclei_A = Nuclei_A - (Nuclei_A / Tau_A) * dt
    Nuclei_B = Nuclei_B + ((Nuclei_A / Tau_A)-(Nuclei_B / Tau_B)) * dt
    exact_A = Nuclei_0 * math.exp(-t[i+1] / Tau_A)
    # Two solutions for Nuclei B for two cases Tau_A = Tau_B and Tau_A =/ Tau_B
    if Tau_A == Tau_B:
            exact_B = (Nuclei_0* t[i+1] * math.exp(-t[i+1] / Tau_B))/ Tau_A
    else:
            exact_B = (Nuclei_0 * Tau_B)*(math.exp(-t[i+1] / Tau_A)-math.exp(-t[i+1] / Tau_B))/(Tau_A -Tau_B)

    # Computing fractional error
    err_A[i+1] = abs((Nuclei_A - exact_A)/exact_A)
    err_B[i+1] = abs((Nuclei_B - exact_B)/exact_B)


pylab.plot(t, err_A,'o-r', label='Error A')
pylab.plot(t, err_B, 'o-b', label='Error B')
pylab.legend()
pylab.xlabel('Time')
pylab.ylabel('Error_A and Error_B')
pylab.title('Error Radioactive Decay')
pylab.grid()
pylab.show()








