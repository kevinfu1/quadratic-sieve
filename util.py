from math import fabs, ceil, sqrt, exp, log, isqrt
import random

def transpose(matrix):
#transpose matrix so columns become rows, makes list comp easier to work with
    return [list(row) for row in zip(*matrix)]

def matvec(A, x):
    """
    Compute the matrix-vector product A*x.
    
    Arguments:
    A -- the input matrix
    x -- the input vector
    
    Returns:
    y -- the product A*x
    """
    n = len(A)
    m = len(x)
    y = [0] * n
    for i in range(n):
        for j in range(m):
            y[i] += A[i][j] * x[j]
    return y

def dot(x, y):
    """
    Compute the dot product of two vectors.
    
    Arguments:
    x -- the first input vector
    y -- the second input vector
    
    Returns:
    z -- the dot product of x and y
    """
    n = len(x)
    z = 0
    for i in range(n):
        z += x[i] * y[i]
    return z

def householder(x):
    """
    Compute the Householder reflection that maps x to a multiple of the first standard basis vector.
    
    Arguments:
    x -- the input vector
    
    Returns:
    H -- the Householder reflection matrix
    """
    n = len(x)
    v = x[:]
    sign = 1 if v[0] >= 0 else -1
    norm = 0
    for i in range(n):
        norm += v[i] * v[i]
    norm = sign * (norm ** 0.5)
    v[0] += norm
    beta = 2 / dot(v, v)
    H = [[beta * v[i] * v[j] for j in range(n)] for i in range(n)]
    for i in range(n):
        H[i][i] -= beta
    return H

def qr(A):
    """
    Compute the QR decomposition of a matrix A using Householder reflections.
    
    Arguments:
    A -- the input matrix
    
    Returns:
    Q -- the orthogonal matrix Q
    R -- the upper triangular matrix R
    """
    n = len(A)
    m = len(A[0])
    Q = [[float(i == j) for j in range(n)] for i in range(n)]
    R = A[:]
    
    for j in range(m):
        v = [R[i][j] for i in range(j, n)]
        H = householder(v)
        Q = matvec(H, Q)
        R = matvec(H, R)
    
    return Q, R

def eig(A, tol=1e-12):
    """
    Compute the eigenvalues and eigenvectors of a symmetric tridiagonal matrix A using the QR algorithm.
    
    Arguments:
    A -- the input symmetric tridiagonal matrix
    tol -- the tolerance for determining convergence (default 1e-12)
    
    Returns:
    eigvals -- a list of the eigenvalues of A
    eigvecs -- a list of the eigenvectors of A
    """
    n = len(A)
    eigvals = [A[i][i] for i in range(n)]
    eigvecs = [[float(i == j) for j in range(n)] for i in range(n)]
    converged = False
    
    while not converged:
        Q, R = qr(A)
        A = matvec(R, Q)
        eigvecs = matvec(Q, eigvecs)
        converged = True
        for i in range(n-1):
            if abs(A[i+1][i]) > tol:
                converged = False
                break
    
    return eigvals, eigvecs

def block_lanczos(A, m, tol=1e-12):
    """
    Compute the left null space of a matrix A using the block Lanczos algorithm.
    
    Arguments:
    A -- the input matrix
    m -- the number of block Lanczos iterations to perform
    tol -- the tolerance for determining linear independence (default 1e-12)
    
    Returns:
    V -- the orthonormal basis for the left null space of A
    """
    n = len(A)
    k = len(A[0])
    V = [[random.uniform(-1, 1) for _ in range(m)] for _ in range(k)]
    Q = qr(V)
    H = [[0] * m for _ in range(m)]
    
    for i in range(m):
        w = matvec(A, Q[:, i])
        for j in range(i):
            H[j][i] = dot(Q[:, j], w)
            w = [w[k] - H[j][i] * Q[k][j] for k in range(n)]
        H[i][i] = sqrt(dot(w, w))
        if H[i][i] < tol:
            break
        Q[:, i+1] = [w[k] / H[i][i] for k in range(n)]
        
    T = [H[i][:i+1] for i in range(i+1)]
    eigvals, eigvecs = eig(T)
    idx = [abs(eigvals[i]) < tol for i in range(len(eigvals))]
    nullspace_basis = [[0] * sum(idx) for _ in range(n)]
    for i in range(n):
        for j in range(sum(idx)):
            nullspace_basis[i][j] = dot(V[k], eigvecs[k][j] * Q[:, k])
    return nullspace_basis