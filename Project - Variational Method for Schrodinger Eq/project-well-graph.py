# project infinite well - graphing
import math
import pylab

##D = pylab.array([[1,2,3,0],[0,1,9,1],[3,1,5,2],[0,2,3,4]])
D = pylab.array([
[  -9.99989758e-01 ,  0.00000000e+00 , -9.62512214e-01,   0.00000000e+00, 6.17024820e-01],
 [  0.00000000e+00 ,  -3.07099182e+00  , 0.00000000e+00 , -3.89594459e+00, 0.00000000e+00],
 [ 2.33425559e-01  , 0.00000000e+00  , 8.98766320e+00  , 0.00000000e+00, -1.30218207e+01],
 [  0.00000000e+00 , 1.58762744e+00  , 0.00000000e+00  , 1.26452091e+01, 0.00000000e+00],
 [  -1.89585931e-02  , 0.00000000e+00 , -5.99363128e+00  , 0.00000000e+00, 2.62659243e+01]])
N = 5 # number of basis states

# defining a function that calculates polynomial basis functions
def basis(y, n):
    z = (y**n)*(y-1)*(y+1)
    return z
# creating x-axis arrays for plotting
X = pylab.zeros(101)
X[0] = -1.0
for i in range(100):
    X[i+1] = X[i]+0.02

for i in range(1,5):
    # creating list for analytic solutions
    Psi = [0.0]* 101
    # creating array for calulated eigenfunctions
    func = pylab.zeros(101)
    for k in range(N):
        func = func + D[k,i-1] * basis(X,k)
    # odd & even analytic solutions
    if i %2 ==0 :
        j = 0
        for x in X:
            Psi[j] = math.sin(i*math.pi*x /2.0)
            j = j+1
    else:
        j = 0
        for x in X:
            Psi[j] = math.cos(i*math.pi*x /2.0)
            j = j+1
    pylab.figure()
    pylab.plot(X, Psi, 'r')
    pylab.plot(X, func, 'b--')
    pylab.xlim([-1.0,1.0])
    pylab.xlabel('x')
    pylab.ylabel('Wavefunction')
    pylab.title('Stationary state for infinite square well, n = %f' %i)
    pylab.show()
            
            































            


