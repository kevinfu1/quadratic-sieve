from math import fabs, ceil, sqrt, exp, log, isqrt
import random
from shanks import tonelli

#O( NLog(Log N))
def sieve_of_eratosthenes(n):
  #""" Returns  a list of primes < n """
  if n <= 2:
    return []
  sieve = [True] * n
  for p in range(3,int(n**0.5)+1,2):
      if sieve[p]:
        for multiple in range(p*p, n, p):
          sieve[multiple] = False
  return  [i for i in range(3,n,2) if sieve[i]]

# Euclid's algorithm
def gcd(a, b):
  if b == 0:
        return a
  elif a >= b:
      return gcd(b,a % b)
  else:
      return gcd(b,a)

def smoothness_bound(N):
    B = exp(sqrt(log(N)*log(log(N)))*0.5)
    return int(B)

# Euler criterion using Legendre symbol
def legendre_symbol(a, p):
    ls = pow(a, (p-1)//2, p)
    return ls if ls == 1 else -1

def factor_base(n, b):
    primes = sieve_of_eratosthenes(b)
    print("primes: ", primes)
    factor_base = []
    for p in primes:
        #print("p = ", p)
        #print("legendre_symbol(n, p): ", legendre_symbol(n, p))
        if legendre_symbol(n, p) == 1:
            factor_base.append(p)
    return [2] + factor_base

def find_smooth(factor_base, N):
    root = isqrt(N) + 1
    print("root = ", root)
# tries to find B-smooth numbers in sieve_seq, using sieving, ie numbers in sieve_seq that have only prime factors in factor_base

    def sieve_prep(N, root):
    # generates a sequence from Y(x) = x^2 - N, starting at x = root 
        sieve_seq = [(root+i)**2 - N for i in range(6 * len(factor_base))]

        return sieve_seq

    sieve_seq = sieve_prep(N, root)
    sieve_list = sieve_seq.copy() # keep a copy of sieve_seq for later
    if factor_base[0] == 2:
        i = 0
        while sieve_list[i] % 2 != 0:
            i += 1
        for j in range(i,len(sieve_list),2): # found the 1st even term, now every other term will also be even
            while sieve_list[j] % 2 == 0: #account for powers of 2
                sieve_list[j] //= 2
    #print("")
    for p in factor_base[1:]: #not including 2
        #print("p = ", p)
        residues = tonelli(N,p) #finds x such that x^2 = N (mod p). There are two start solutions
        #print("residues: ", residues)
        for r in residues:
            for i in range((r-root) % p, len(sieve_list), p): # Now every pth term will also be divisible
                while sieve_list[i] % p == 0: #account for prime powers
                    sieve_list[i] //= p
                    
    #xlist = [] #original x terms
    print("sieve_seq: ", sieve_seq)
    print("sieve_list: ", sieve_list)
    smooth_nums = []
    #indices = [] # index of discovery
    
    for i in range(len(sieve_list)):
        if len(smooth_nums) >= len(factor_base): # found enough B-smooth numbers
            break
        if sieve_list[i] == 1: # found B-smooth number
            smooth_nums.append(sieve_seq[i])
           # xlist.append(i+root)
            #indices.append(i)

    return smooth_nums
  

def build_matrix(smooth_nums,factor_base):
# generates exponent vectors mod 2 from previously obtained smooth numbers, then builds matrix

    def factor(n,factor_base):#trial division from factor base
        factors = []
        if n < 0:
            factors.append(-1)
        for p in factor_base:
            if p == -1:
                pass
            else:
                while n % p == 0:
                    factors.append(p)
                    n //= p
        return factors

    M = []
    factor_base.insert(0,-1)
    for n in smooth_nums:
        exp_vector = [0]*(len(factor_base))
        n_factors = factor(n,factor_base)
        #print(n,n_factors)
        for i in range(len(factor_base)):
            if factor_base[i] in n_factors:
                exp_vector[i] = (exp_vector[i] + n_factors.count(factor_base[i])) % 2

        #print(n_factors, exp_vector)
        if 1 not in exp_vector: #search for squares
            return True, n
        else:
            pass
        
        M.append(exp_vector)  
    #print("Matrix built:")
    #mprint(M)
    return(False, transpose(M))

    
def transpose(matrix):
#transpose matrix so columns become rows, makes list comp easier to work with
    new_matrix = []
    for i in range(len(matrix[0])):
        new_row = []
        for row in matrix:
            new_row.append(row[i])
        new_matrix.append(new_row)
    return(new_matrix)

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



def quad_sieve(n):
  B = smoothness_bound(n)
  print("B = ", B)
  factor_b = factor_base(n, B)
  print("factor_base: ", factor_b)
  smooth_nums = find_smooth(factor_b, n)
  print("smooth_nums: ", smooth_nums)

  print("Building matrix...")


if __name__ == "__main__":
  quad_sieve(227179)