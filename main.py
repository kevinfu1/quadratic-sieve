from math import fabs, ceil, sqrt, exp, log, isqrt
import random, time
from shanks import tonelli
from helper import gauss_elim, solve_row, solve, transpose, gcd


def sieve_of_eratosthenes(n):
    """
    This function returns a list of primes less than n using the Sieve of Eratosthenes algorithm.

    :param n: The upper bound for the list of primes.
    :return: A list of primes less than n.

    Time complexity: O(n * log(log(n)))
    """

    if n <= 2:  # If n is less than or equal to 2, there are no primes less than n
        return []
    sieve = [True] * n  # Creating a boolean list of size n, initialized to True
    for p in range(3, int(n**0.5) + 1, 2):  # Iterating through odd numbers up to the square root of n
        if sieve[p]:  # If the current number is prime
            for multiple in range(p * p, n, p):  # Marking the multiples of the prime number as composite
                sieve[multiple] = False
    # Returning the list of prime numbers, skipping the even numbers
    return [i for i in range(3, n, 2) if sieve[i]]


def smoothness_bound(N):
    """
    This function computes the smoothness bound for the quadratic sieve algorithm.

    :param N: The number to be factored.
    :return: The smoothness bound.

    Time complexity: O(1)
    """

    B = exp(sqrt(log(N) * log(log(N))) * 0.5)  # Computing the smoothness bound using the formula B = exp(sqrt(log(N) * log(log(N))) * 0.5)
    return int(B)  # Returning the smoothness bound as an integer

def legendre_symbol(a, p):
    """
    This function uses the Euler criterion to compute the Legendre symbol (a/p).

    :param a: The number whose symbol is to be computed.
    :param p: The prime number.
    :return: The Legendre symbol (a/p).

    Time complexity: O(log(p)), where p is the prime number
    """

    ls = pow(a, (p - 1) // 2, p)  # Computing the Legendre symbol using the Euler criterion
    return ls if ls == 1 else -1  # Returning the symbol, which is either 1 or -1 depending on the result of the computation


# Time complexity: O(n * log(log(n))) due to sieve_of_eratosthenes call
def factor_base(n, b):
    """
    This function computes the factor base for the quadratic sieve algorithm.

    :param n: The number to be factored.
    :param b: The size of the factor base.
    :return: A list containing the primes in the factor base.
    """

    primes = sieve_of_eratosthenes(b)  # Generating a list of primes up to b using the sieve of Eratosthenes
    factor_base = [2]  # Adding -1 and 2 to the factor base
    for p in primes:
        if legendre_symbol(n, p) == 1:  # Checking if n is a quadratic residue modulo p
            factor_base.append(p)  # Adding p to the factor base if it meets the criteria

    return factor_base


def find_smooth(factor_base, N, I, root):
    """
    This function finds B-smooth numbers for the quadratic sieve algorithm.

    :param factor_base: A list of primes to be used as the factor base.
    :param N: The number to be factored.
    :param I: The size of the sieve interval.
    :return: A tuple containing two lists. The first list contains the B-smooth numbers found, and the second
             list contains the corresponding x values that were used to find those B-smooth numbers.
    """


    def sieve_prep(N, root, sieve_int):
        """
        This function generates a sequence from Y(x) = x^2 - N, starting at x = root, which will be used for the sieve.

        :param N: The number to be factored.
        :param root: The integer square root of N plus 1.
        :param sieve_int: The size of the sieve interval.
        :return: A list containing the sieve sequence and the last number checked in this sieve sequence


        """

        sieve_seq = [x**2 - N for x in range(root,root+sieve_int)] # Generating the sieve sequence
        return sieve_seq, root+sieve_int

    # Prepare the sieve sequence
    sieve_seq, new_root = sieve_prep(N, root, I) # Generating the sieve sequence
    sieve_list = sieve_seq.copy() # Making a copy of sieve_seq for later

    # Dividing out powers of 2 from the sieve sequence
    if factor_base[0] == 2:
        i = 0
        while sieve_list[i] % 2 != 0:
            i += 1
        for j in range(i,len(sieve_list),2): # found the 1st even term, now every other term will also be even
            while sieve_list[j] % 2 == 0: #account for powers of 2
                sieve_list[j] //= 2
   
    # Dividing out powers of primes from the sieve sequence
    for p in factor_base[1:]:  # Looping through the factor base, not including 2
        residues = tonelli(N,p)  # Finding the square roots of N modulo p
        for r in residues:
            for i in range((r-root) % p, len(sieve_list), p):
                while sieve_list[i] % p == 0:
                    sieve_list[i] //= p
                    
    # Finding the B-smooth numbers in the sieve sequence
    xlist = []  # Creating an empty list to hold the corresponding original x values
    smooth_nums = [] # Creating an empty list to hold the B-smooth numbers found
    
    for i in range(len(sieve_list)):
        if len(smooth_nums) >= len(factor_base) + 1:  # If enough B-smooth numbers have been found, break out of the loop
            break
        if sieve_list[i] == 1:  # If the sieve value is 1, the corresponding number is B-smooth
            smooth_nums.append(sieve_seq[i])  # Adding the B-smooth number to the list
            xlist.append(i+root)  # Adding the corresponding x value to the list

    return smooth_nums, xlist, new_root
  
def factor(n, factor_base):
    """
    This function factorizes a number using the given factor base.

    :param n: The number to be factorized.
    :param factor_base: The factor base to be used.
    :return: A list of prime factors of n.

    Time complexity: The time complexity of the function is O(k log n), n is the number to be factorized and k is the size of the factor base
    """

    factors = []
    if n < 0:
        factors.append(-1)
        n = -n
    for p in factor_base:  # Checking if n is divisible by each prime in the factor base
        if p != -1 and n % p == 0:  # If p is a factor of n
            while n % p == 0:  # Divide n by p until it is no longer divisible by p
                factors.append(p)  # Add p to the list of factors
                n //= p
    if n > 1:  # If n is still greater than 1, it must be prime and is added to the list of factors
        factors.append(n)
    return factors


# Builds a matrix for the quadratic sieve algorithm
def build_matrix(smooth_nums, factor_base):
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

# Handles the case when a perfect square is found
def handleSquare(M, n, smooth_nums, xlist):
    factor_1 = None
    factor_2 = None
    x = smooth_nums.index(M)
    print(x)
    print(xlist)
    factor_1 = gcd(xlist[x] + sqrt(M), n)
    factor_2 = n // factor_1
    return factor_1, factor_2

def quad_sieve(n, I):
  B = smoothness_bound(n)
  print("B = ", B)
  factor_b = factor_base(n, B)
  #print("factor_base: ", factor_b)
  smooth_nums = []
  xlist = []
  root = isqrt(n) + 1
  while len(smooth_nums) < len(factor_b):
    smooth, x, root = find_smooth(factor_b, n, I, root)
    smooth_nums += smooth
    xlist += x

  #smooth_nums, xlist = find_smooth(factor_b, n, I)
  print("smooth_nums: ", smooth_nums)


  if len(smooth_nums) < len(factor_b):
    return("Not enough smooth numbers. Increase the sieve interval or size of the factor base.")

  print("Building matrix...")
  is_square, t_matrix = build_matrix(smooth_nums,factor_b)
  print("Matrix built: ", t_matrix)
  
  if is_square:
        factor_1, factor_2 = handleSquare(t_matrix, n, smooth_nums, xlist)
        print('square found')
        return factor_1, factor_2
  else:
    sol_rows, marks, M = gauss_elim(t_matrix)
    print("M: ", M)
    print("sol_rows: ", sol_rows)
    solution_vec = solve_row(sol_rows, M, marks, 0)
    print("solution_vec: ", solution_vec)
    factor = solve(solution_vec, smooth_nums, xlist, n)
    print("factor: ", factor)
    for K in range(1,len(sol_rows)):
      if (factor == 1 or factor == n):
          print("Didn't work. Trying different solution vector...")
          solution_vec = solve_row(sol_rows,M,marks,K)
          factor = solve(solution_vec,smooth_nums,xlist,n)
      else:
          print("Found factors!")
          return factor, int(n/factor)
  
  return("Didn't find any nontrivial factors!")


if __name__ == "__main__":
    start_time = time.time()
    #print(quad_sieve(16921456439215439701,10000))
  
#   print(quad_sieve(46839566299936919234246726809, pow(10, 4)))

    #print(quad_sieve(6172835808641975203638304919691358469663, pow(10, 4)))

    #print(quad_sieve(3744843080529615909019181510330554205500926021947, 10000))
    #print(quad_sieve(125513, 10000))

    print("--- %s seconds ---" % (time.time() - start_time))