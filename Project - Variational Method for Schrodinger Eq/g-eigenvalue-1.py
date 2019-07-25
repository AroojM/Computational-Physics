# generalized eigenvalue
# example - choleski + jacobi
# jacobi as function
import math
import pylab
# generalized eignevalue problem H $\alpha$ = E S $\alpha$

H = pylab.array([[1./3., -1./3, 0.],[-1./3., 4./3., -1.],[0., -1., 2.]])
S = pylab.array([[1., 0., 0.],[0., 1., 0.],[0., 0., 2.]])
#--------------------------------------------------------------
# Defining Functions ------------------------------------------
#--------------------------------------------------------------
# defining a function that carries out Cholesky decomposition
# & returns the lower triangular matrix
def cholesky(M):
    n = len(M)
    L = pylab.zeros((n,n), dtype = float)
    for j in range(n):
        L[j,j] = math.sqrt(M[j,j]- pylab.dot(M[j, 0:j],M[j, 0:j]))
    for i in range(j+1, n):
        L[i,j] = (M[i,j] - pylab.dot(M[i,0:j],M[j,0:j]))/M[j,j]
    return L

# defining a function that inverts a lower triangular matrix
# by forward substitution
def inverse(M):
    n = len(M)
    for j in range(n-1):
        M[j,j] = 1.0/M[j,j]
        for i in range(j+1,n):
            M[i,j] = - pylab.dot(M[i,j:i],M[j:i,j])/ M[i,i]
    M[n-1, n-1] = 1.0/M[n-1,n-1]
    
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
#--------------------------------------------------------------
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
                phi = aDiff/(2.0*a[p,q])
                t = 1.0/(abs(phi) + math.sqrt(phi**2 + 1.0))
                if phi < 0.0:
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
            for r in range(q+1,n1):  # For r > q
                temp = a[p,r]
                a[p,r] = temp - s*(a[q,r] + tau*temp)
                a[q,r] = a[q,r] + s*(temp - tau*a[q,r])
            for r in range(n1):      # Updating transformation matrix
                temp = v[r,p]
                v[r,p] = c*temp - s*v[r,q]
                v[r,q] = s*temp + c*v[r,q]
        if g == maxRot-1 and aMax > tol : print 'Jacobi method needs more rotation matrices to reach the desired tolerance level'

    Eval = pylab.diagonal(a) # eigenvalues 
    Evec = v           # eigenvectors 
    return Eval, Evec
    # end - Jacobi diagonalization

#--------------------------------------------------------------
#-----------------------------------------------------------

# cholesky decomposition of S
L = cholesky(S)

# inverting L
inverse(L)

# hamiltonian matrix when equation reduced to standard eigenvalue problem
Hbar = pylab.dot(L,pylab.inner(H,L))

# transpose of L inverse
C = pylab.transpose(L) 

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

print 'Eigenvalues of Hbar'
print Eval
print 'Eigenvectors of Hbar'
print Evec
# Eigenvectors - 'alpha' of H
D = pylab.dot(C, Evec)
print 'Eigenvalues of H'
print Eval
print 'Eigenvectors of H'
print D












