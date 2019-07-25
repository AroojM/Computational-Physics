
# Exercise 1.1

#Exercise 1.1 Produce a program that is similar to plotcos.py but do not
#use the cos function from the math module (or any other module). Instead,
#compute it using ordinary multiplication and additions.



import pylab

steps = 200
dtheta = 0.1

# define function fact(x) that computes x!
def fact(x):
    z = 1
    y = x+1
    for k in range(1,y):
        z = k * z
    return z

# create lists for the values of theta and cosine
theta = [0.0] * steps
cosine = [0.0] * steps

# nested loops
# j loop fills in values of theta and cosine
# i loop computes values of cosine by Taylor series
for j in range(steps):
    theta[j] = j * dtheta
    cosine[j] = 0.0
    for i in range (0,40):
        cosine[j]= ((((-1)**i)*((theta[j]**(2*i))))/fact(2*i)) + cosine[j]
    

pylab.plot(theta, cosine, 's-b', label='cos')
pylab.xlabel('theta(rad)')
pylab.ylabel('cos(theta)')
pylab.legend()
pylab.title('Plot of cosine function (from Taylor series)')
pylab.grid()
pylab.show()
