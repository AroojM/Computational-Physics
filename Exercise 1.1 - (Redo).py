# Exercise 1.1 (2nd method)

import pylab

steps = 200
dtheta = 0.1

# create lists for the values of theta, cosine and angle
theta = [0.0] * steps
cosine = [0.0] * steps
angle = [0.0] * steps

# define COS function using first four terms of Maclaurian series
def COS(x):
    y = 1.0 - ((x**2)/2.0) + ((x**4)/24.0) - ((x**6)/720.0)
    return y

# j loop fills in values of theta, cosine and angle
for j in range(steps):
    theta[j] = j * dtheta
    # angle maps values of theta greater than 2pi to (0,2pi) range
    angle[j] = theta[j] % 6.28318530718
    # using chained conitional to compute COS for different ranges of angle
    if angle[j] < 1.57079632679:
        cosine[j] = COS(angle[j])
    elif 1.57079632679 < angle[j] < 4.71238898038:
        cosine[j] = -COS(3.14159265359 - angle[j])
    else:
        cosine[j] = COS(6.28318530718-angle[j])

# plot COS vs theta
pylab.plot(theta, cosine, 'o-g', label='cos')
pylab.xlabel('theta(rad)')
pylab.ylabel('cos(theta)')
pylab.legend()
pylab.title('Plot of cosine function')
pylab.grid()
pylab.show()
