from math import sqrt, exp, log, isqrt
import time
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

    Time complexity: O(b log^2(b))
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
    :param root: The square root of the number N being factored, rounded up to the nearest integer.
    :return: A tuple containing two lists. The first list contains the B-smooth numbers found, and the second
             list contains the corresponding x values that were used to find those B-smooth numbers.

    Time Complexity: O(I k log^2 n), where I is the number of smooth numbers to be found, k is the size of the factor base, and n is the number being factored
    """


    def sieve_prep(N, root, sieve_int):
        """
        This function generates a sequence from Y(x) = x^2 - N, starting at x = root, which will be used for the sieve.

        :param N: The number to be factored.
        :param root: The integer square root of N plus 1.
        :param sieve_int: The size of the sieve interval.
        :return: A list containing the sieve sequence and the last number checked in this sieve sequence

        Time Complexity: O(sieve_int) 
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

    Time complexity: O(k log n), where k is the size of the factor base
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


def build_matrix(smooth_nums, factor_base):
    """
    This function builds a matrix for the quadratic sieve algorithm.

    :param smooth_nums: A list of smooth numbers.
    :param factor_base: The factor base to be used.
    :return: A tuple containing a Boolean flag indicating if a square has been found,
             and the transpose of the matrix.

    Time complexity: O(N * B^2), where N is the number of smooth numbers and B is the size of the factor base.
    """

    M = []  # Initialize an empty matrix M
    square_found = False  # Initialize the square-found flag to False
    for n in smooth_nums:
        exp_vector = [0] * len(factor_base)  # Initialize an exponent vector of zeros for this smooth number
        n_factors = factor(n, factor_base)  # Compute the factorization of n with respect to the factor base
        for i, p in enumerate(factor_base):
            if p in n_factors:
                exp_vector[i] = (exp_vector[i] + n_factors.count(p)) % 2  # Set the corresponding entry in the exponent vector to 1 if the exponent is odd
        M.append(exp_vector)  # Add the exponent vector as a row to the matrix
        if not square_found and all(e == 0 for e in exp_vector):
            square_found = True  # If the row contains all zeros, a square has been found in the matrix
            return square_found, n  # Return the square-found flag and the smooth number that caused it
    
    return square_found, transpose(M)  # Return the square-found flag and the transpose of the matrix

def handleSquare(M, n, smooth_nums, xlist):
    """
    This function handles the case when a perfect square is found in the matrix.

    :param M: The smooth number that is a perfect square.
    :param n: The number to be factored.
    :param smooth_nums: A list of smooth numbers.
    :param xlist: A list of the x values corresponding to the smooth numbers.
    :return: A tuple containing the two factors of n.

    Time complexity: O(1).
    """

    factor_1 = None  # Initialize factor 1 to None
    factor_2 = None  # Initialize factor 2 to None
    x = smooth_nums.index(M)  # Get the index of the smooth number that is a perfect square
    factor_1 = gcd(xlist[x] + sqrt(M), n)  # Compute one factor of n using the quadratic sieve algorithm
    factor_2 = n // factor_1  # Compute the other factor of n by dividing n by the first factor
    
    return factor_1, factor_2  # Return the two factors of n as a tuple

def generate_xs(smooth_nums, factor_base, N):
    """
    This function generates a list of x values for the smooth numbers in the matrix.

    :param smooth_nums: A list of smooth numbers.
    :param factor_base: The factor base to be used.
    :param N: The number to be factored.
    :return: A list of x values corresponding to the smooth numbers.

    Time complexity: O(k^2 * log(N)), where k is the size of the factor base.
    """
    root = isqrt(N) + 1  # Compute the square root of N and add 1 to get the starting value of x
    x_list = []  # Initialize the list of x values to an empty list
    for n in smooth_nums:  # Iterate through the list of smooth numbers
        x = root  # Initialize x to the starting value
        for p in factor(n, factor_base):  # Iterate through the prime factors of the smooth number
            if p != -1:
                # While x is divisible by p and (x^2 - N) / n is divisible by p, divide x by p
                while x % p == 0 and (x**2 - N) // n % p == 0:
                    x //= p
        x_list.append(x)  # Add the final value of x to the list
    return x_list  # Return the list of x values


def quad_sieve(n, I):
    """
    Factorizes n using the quadratic sieve algorithm.

    :param n (int): the integer to be factored.
    :param I (int): the size of the sieve interval.
    :return: A tuple of nontrivial factors of n if found, or a message indicating that no nontrivial factors were found.

    Time Complexity: O(B^3), where B is the smoothness bound
    """

    B = smoothness_bound(n)  # Calculate a bound B on the size of smooth numbers to search for. O(1)
    # print("B = ", B)
    factor_b = factor_base(n, B)  # Calculate the factor base of n up to the bound B. O(B * log(log(B)))

    smooth_nums = []
    xlist = []
    root = isqrt(n) + 1  # Find the square root of n, rounded up to the nearest integer. O(log(n))
    while len(smooth_nums) < len(factor_b):
        # Find smooth numbers in the interval (root, root + I) using the factor base. O(B * sqrt(n) * log(n))
        smooth, x, root = find_smooth(factor_b, n, I, root)
        smooth_nums += smooth
        xlist += x

    # print("smooth_nums: ", smooth_nums)

    if len(smooth_nums) < len(factor_b):
        # If there are not enough smooth numbers, return an error message. O(1)
        return "Not enough smooth numbers. Increase the sieve interval or size of the factor base."

    # print("Building matrix...")
    is_square, t_matrix = build_matrix(smooth_nums, factor_b)  # Build the matrix of relations. O(B^2 * log(n))

    if is_square:
        # If the matrix is square, handle the case where one of the factors is a perfect square. O(B * log(n))
        xs = generate_xs(smooth_nums, factor_b, n)
        factor_1, factor_2 = handleSquare(t_matrix, n, smooth_nums, xs)
    else:
        # Otherwise, perform Gaussian elimination on the matrix to find a solution vector. O(B^3)
        sol_rows, marks, M = gauss_elim(t_matrix)
        # print("M: ", M)
        # print("sol_rows: ", sol_rows)
        # Solve the system of equations using the solution vector. O(B * log(n))
        solution_vec = solve_row(sol_rows, M, marks, 0)
        # print("solution_vec: ", solution_vec)
        factor = solve(solution_vec, smooth_nums, xlist, n)
        # print("factor: ", factor)
        for K in range(1, len(sol_rows)):
            if factor in (1, n):
                # If the factor is 1 or n, try a different solution vector. O(B^2)
                # print("Didn't work. Trying different solution vector...")
                solution_vec = solve_row(sol_rows, M, marks, K)
                factor = solve(solution_vec, smooth_nums, xlist, n)
            else:
                # If a nontrivial factor is found, return it. O(1)
                # print("Found factors!")
                return factor, n // factor

    # If no nontrivial factors are found, return an error message. O(1)
    return "Didn't find any nontrivial factors!"

if __name__ == "__main__":
    start_time = time.time()

    # print(quad_sieve(16921456439215439701,pow(10, 4)))
    print(quad_sieve(46839566299936919234246726809, pow(10, 4)))
    # print(quad_sieve(6172835808641975203638304919691358469663, pow(10, 4)))
    # print(quad_sieve(3744843080529615909019181510330554205500926021947, pow(10, 4)))

    print("--- %s seconds ---" % (time.time() - start_time))