# project hydrogen atom


import math
import pylab

# creating matrices
gamma = [13.00773, 1.962079, 0.444529, 0.1219492]

S = pylab.zeros((4,4), dtype = float)
H = pylab.zeros((4,4), dtype = float)

for m in range(4):
    for n in range(4):
        S[m,n] = (math.pi / (gamma[m]+gamma[n]))**(3.0/2.0)
        H[m,n] = ((3*gamma[m]*gamma[n]*(math.pi)**(3.0/2.0))/(gamma[m]+gamma[n])**(5.0/2.0))-(2*math.pi/(gamma[m]+gamma[n]))
#------------------------------------------------------------------------------
# Defining Functions ----------------------------------------------------------
#------------------------------------------------------------------------------

# defining a function that inverts a matrix
#  might be unstable if division by zero - in such case use pivoting method
def inverse(A):
    def LUdecomp(A):
        # use Crout's algorithm to perform LU decomposition of A
        n = len(A)
        L = pylab.zeros(A.shape)
        U = pylab.zeros(A.shape)
        for i in range(n): L[i,i] = 1.0
        for j in range(n):
            for i in range(j+1):
                U[i,j] = A[i,j]
                for k in range(i):
                    U[i,j] -= L[i,k] * U[k,j]
            for i in range(j+1,n):
                L[i,j] = A[i,j]
                for k in range(j):
                    L[i,j] -= L[i,k] * U[k,j]
                L[i,j] /= U[j,j]
        return L, U

    def solve(A,b):
        # solves the linear system A.x = b for x
        n = len(A)
        L, U = LUdecomp(A)
        x = pylab.zeros(b.shape)
        y = pylab.zeros(b.shape)
        # forward substitute to solve equation L.y = b for y
        for i in range(n):
            y[i] = b[i]
            for j in range(i):
                y[i] -= L[i,j] * y[j]
            y[i] /= L[i,i]
        # back substitute to solve equation U.x = y for x
        for i in reversed(range(n)):
            x[i] = y[i]
            for j in range(i+1,n):
                x[i] -= U[i,j] * x[j]
            x[i] /= U[i,i]
        return x
    # when b is a matrix, solve(A,b) gives inverse of A
    B = pylab.eye(len(A))
    return solve(A,B)
#-----------------------------------------------------------------------------   
# defining swap functions
# swap1 - swaps entries in vector
def swap1(V, i, j):
    temp = pylab.copy(V[i])
    V[i] = V[j]
    V[j] = temp
# swap2 - swaps columns in martix
def swap2 (M, i, j):
    temp = pylab.copy(M[:,j])
    M[:,j] = M[:,i]
    M[:,i] = temp
#-----------------------------------------------------------------------------
# defining a function that carries out jacobi transformation
def jacobi (a, n):    
    # tolerance
    tol = 1.0e-9
    # transformation matrix initilazied to identity matrix
    v = pylab.identity(n)*1.0
    # maximum number of Jacobi rotations
    maxRot = 5*(n**2)

    # Jacobi method - eigenvalues & eigenvectors 
    g = 0
    for g in range(maxRot):
        # searching for largest element
        aMax = 0.0
        for r in range(n-1):
            for u in range(r+1,n):
                if abs(a[r,u]) >= aMax:
                    aMax = abs(a[r,u])
                    p = r; q = u
        # start - Jacobi diagonalization 
        if aMax > tol:
            aDiff = a[q,q] - a[p,p]
            if abs(a[p,q]) < abs(aDiff)*1.0e-40:
                t = a[p,q]/aDiff
            else:
                theta = aDiff/(2.0*a[p,q])
                t = 1.0/(abs(theta) + math.sqrt(theta**2 + 1.0))
                if theta < 0.0:
                    t = -t
            c = 1.0/math.sqrt(t**2 + 1.0); s = t*c
            tau = s/(1.0 + c)
            temp = a[p,q]
            a[p,q] = 0.0
            a[p,p] = a[p,p] - t*temp
            a[q,q] = a[q,q] + t*temp
            for r in range(p):      # For r < p
                temp = a[r,p]
                a[r,p] = temp - s*(a[r,q] + tau*temp)
                a[r,q] = a[r,q] + s*(temp - tau*a[r,q])
            for r in range(p+1,q):  # For p < r < q
                temp = a[p,r]
                a[p,r] = temp - s*(a[r,q] + tau*a[p,r])
                a[r,q] = a[r,q] + s*(temp - tau*a[r,q])
            for r in range(q+1,n):  # For r > q
                temp = a[p,r]
                a[p,r] = temp - s*(a[q,r] + tau*temp)
                a[q,r] = a[q,r] + s*(temp - tau*a[q,r])
            for r in range(n):      # Updating transformation matrix
                temp = v[r,p]
                v[r,p] = c*temp - s*v[r,q]
                v[r,q] = s*temp + c*v[r,q]
        if g == maxRot-1 and aMax > tol : print 'Jacobi method needs more rotation matrices to reach the desired tolerance level'

    Eval = pylab.diagonal(a) # eigenvalues 
    Evec = v           # eigenvectors 
    return Eval, Evec
#-------------------------------------------------------------------------------
# End - defining functions -----------------------------------------------------
#-------------------------------------------------------------------------------

h = len(S)
s, U = jacobi(S, h)
iss = pylab.zeros((h,h), dtype = float) # inverse square-root of s
for i in range (h):
    for j in range(h):
        if j == i :
            iss[i,j] = 1.0 / (math.sqrt(s[i]))
V = pylab.dot(U, iss)

# hamiltonian matrix when equation reduced to standard eigenvalue problem
Hbar = pylab.dot(pylab.transpose(V),pylab.dot(H,V))
n1 = len(Hbar)
Eval, Evec = jacobi(Hbar, n1)

# sorting output from Jacobi diagonalization
n2 = len(Eval)
for i in range(n2-1):
    index = i
    val = Eval[i]
    for j in range(i+1, n2):
        if Eval[j] < val:
            index = j
            val = Eval[j]
    if index != i:
        swap1(Eval, i, index)
        swap2(Evec, i, index)

# Eigenvalues of H and Hbar - Eval
# Eigenvectors - 'beta' of Hbar - Evec

inverse(V)

# Eigenvectors - 'alpha' of H
D = pylab.dot(V, Evec)

print 'Ground state energy of hydrogen atom =', Eval[0],'Hartees = ',Eval[0]*27.212, 'eV'




