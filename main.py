from math import fabs, ceil, sqrt, exp, log, isqrt
import random
from shanks import tonelli
from util import block_lanczos, transpose

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

#changed factor base to add -1 as possibility as well incase number is negative
def factor_base(n, b):
    primes = sieve_of_eratosthenes(b)
    factor_base = [-1, 2]  # add -1 and 2 to the factor base
    for p in primes:
        if legendre_symbol(n, p) == 1:
            factor_base.append(p)
    return factor_base

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
  
def factor(n, factor_base):
    factors = []
    if n < 0:
        factors.append(-1)
        n = -n
    for p in factor_base:
        if p != -1 and n % p == 0:
            while n % p == 0:
                factors.append(p)
                n //= p
    if n > 1:
        factors.append(n)
    return factors

def build_matrix(smooth_nums,factor_base):
# generates exponent vectors mod 2 from previously obtained smooth numbers, then builds matrix

    M = []
    square_found = False
    for n in smooth_nums:
        exp_vector = [0] * len(factor_base)
        n_factors = factor(n, factor_base)
        for i, p in enumerate(factor_base):
            if p in n_factors:
                exp_vector[i] = (exp_vector[i] + n_factors.count(p)) % 2
        M.append(exp_vector)
        if not square_found and all(e == 0 for e in exp_vector):
            square_found = True
            return square_found, n
    return square_found, transpose(M)

def handleSquare(M, n, smooth_nums, xlist):
    factor_1 = None
    factor_2 = None
    x = smooth_nums.index(M)
    factor_1 = gcd(xlist[x]+sqrt(M),n)
    factor_2 = n // factor_1   
    return factor_1, factor_2

def generate_xs(smooth_nums, factor_base, N):
    root = isqrt(N) + 1
    x_list = []
    for n in smooth_nums:
        x = root
        for p in factor(n, factor_base):
            if p != -1:
                while x % p == 0 and (x**2 - N) // n % p == 0:
                    x //= p
        x_list.append(x)
    return x_list

def solve(t_matrix, is_square, n, smooth_nums, factor_b):
    if is_square:
        xs = generate_xs(smooth_nums, factor_b, n)
        factor_1, factor_2 = handleSquare(t_matrix, n, smooth_nums, xs)
    else:
        #solve matrix, how many iterations?
        ln = block_lanczos(t_matrix, 10)
    return factor_1, factor_2

def quad_sieve(n):
  B = smoothness_bound(n)
  print("B = ", B)
  factor_b = factor_base(n, B)
  print("factor_base: ", factor_b)
  smooth_nums = find_smooth(factor_b, n)
  print("smooth_nums: ", smooth_nums)

  print("Building matrix...")
  is_square, t_matrix = build_matrix(smooth_nums,factor_b)

  factor_1, factor_2 = solve(t_matrix, is_square, n, smooth_nums, factor_b)
  return factor_1, factor_2


if __name__ == "__main__":
  quad_sieve(227179)