from math import fabs, ceil, sqrt, exp, log, isqrt
import random
from shanks import STonelli

def sieve_of_eratosthenes(n):
  #""" Returns  a list of primes < n """
  if n <= 2:
    return []
  sieve = [True] * n
  for p in range(3,int(n**0.5)+1,2):
      if sieve[p]:
        for multiple in range(p*p, n, p):
          sieve[multiple] = False
  return [2] + [i for i in range(3,n,2) if sieve[i]]

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
    #I = F**3
    #print(F,I)
    return int(B)

# Euler criterion using Legendre symbol
def legendre_symbol(a, p):
    ls = pow(a, (p-1)//2, p)
    return ls if ls == 1 else -1

def factor_base(n):
    primes = sieve_of_eratosthenes(n)
    factor_base = []
    for p in primes:
        if legendre_symbol(n, p) == 1:
            factor_base.append(p)
    return factor_base

def find_smooth(factor_base, N):
# tries to find B-smooth numbers in sieve_seq, using sieving

    def sieve_prep(N):
    # generates a sequence from Y(x) = x^2 - N, starting at x = root 
        x = isqrt(N) + 1
        sieve_seq = [(x+i)**2 - N for i in range(len(factor_base))]

        return sieve_seq

    sieve_seq = sieve_prep(N)
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
        residues = STonelli(N,p) #finds x such that x^2 = n (mod p). There are two start solutions
        
        for r in residues:
            for i in range((r-root) % p, len(sieve_list), p): # Now every pth term will also be divisible
                while sieve_list[i] % p == 0: #account for prime powers
                    sieve_list[i] //= p
    xlist = [] #original x terms
    smooth_nums = []
    indices = [] # index of discovery
    
    for i in range(len(sieve_list)):
        if len(smooth_nums) >= len(factor_base)+T: #probability of no solutions is 2^-T
            break
        if sieve_list[i] == 1: # found B-smooth number
            smooth_nums.append(sieve_seq[i])
            xlist.append(i+root)
            indices.append(i)

    return(smooth_nums,xlist,indices)
  
# def B_smooth_values(n, factor_base):
#   B_smooth_values = []
#   x = isqrt(n) + 1
#   for i in range(x, len(factor_base)+1):
#     # if gcd(i, n) == 1:
#     #   smooth_values.append(i)

def quad_sieve(n):
  B = smoothness_bound(n)
  factor_base = factor_base(B)



# def isqrt(n): # Newton's method, returns exact int for large squares
#     x = n
#     y = (x + 1) // 2
#     while y < x:
#         x = y
#         y = (x + n // x) // 2
#     return x

# def sieve(n, b):
#     primes = [i for i in range(2, b+1) if is_prime(i)]
#     factor_base = []
#     for p in primes:
#         if legendre_symbol(n, p) == 1:
#             e = 0
#             while n % p == 0:
#                 n //= p
#                 e += 1
#             factor_base.append((p, e))
#     if n > 1:
#         return None
#     return factor_base

# def quad_residue(a, n, p):
#     x = pow(a, (p+1)//4, p)
#     return x if pow(x, 2, p) == a else n-x



# def smooth_sieve(n, b, factor_base):
#     sieve = []
#     for i in range(1, b+1):
#         x = quad_residue(i, n, factor_base[0])
#         factors = []
#         for p in factor_base:
#             e = 0
#             while x % p == 0:
#                 x //= p
#                 e += 1
#             factors.append(e)
#         sieve.append((i, factors))
#     return sieve

# def matrix_rank(matrix):
#     # Implement Gaussian elimination to determine the rank of the matrix
#     m, n = len(matrix), len(matrix[0])
#     rank = 0
#     for j in range(n):
#         for i in range(rank, m):
#             if matrix[i][j] != 0:
#                 for k in range(j, n):
#                     matrix[i][k], matrix[rank][k] = matrix[rank][k], matrix[i][k]
#                 for u in range(rank+1, m):
#                     if matrix[u][j] != 0:
#                         c = matrix[u][j] * pow(matrix[rank][j], -1, n)
#                         for v in range(j, n):
#                             matrix[u][v] = (matrix[u][v] - c * matrix[rank][v]) % n
#                 rank += 1
               
