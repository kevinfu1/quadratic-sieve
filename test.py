from math import isqrt
from random import randrange

# implement quadratic sieve algorithm
def quadratic_sieve(n):
    # find smoothness bound
    B = smoothness_bound(n)
    # find factor base
    factor_base = sieve(n, B)
    # find smooth numbers
    smooth_numbers = smooth_sieve(n, B, factor_base)
    # find linear relations
    matrix = []
    for i in range(len(smooth_numbers)):
        matrix.append(smooth_numbers[i][1])
    # find nullspace
    nullspace = find_nullspace(matrix)
    # find smooth relations
    smooth_relations = []
    for v in nullspace:
        x = 1
        for i in range(len(v)):
            x *= pow(factor_base[i], v[i], n)
        x = x % n
        if x > n/2:
            x = n - x
        if x != 1 and x != n-1:
            smooth_relations.append(x)
    # find factors
    factors = []
    for x in smooth_relations:
        factors.append(gcd(x-1, n))
        factors.append(gcd(x+1, n))
    return factors

def find_nullspace(matrix):
    # find nullspace of matrix
    # using gaussian elimination
    # https://rosettacode.org/wiki/Reduced_row_echelon_form#Python
    lead = 0
    rowCount = len(matrix)
    columnCount = len(matrix[0])
    for r in range(rowCount):
        if lead >= columnCount:
            return
        i = r
        while matrix[i][lead] == 0:
            i += 1
            if i == rowCount:
                i = r
                lead += 1
                if columnCount == lead:
                    return
        matrix[i], matrix[r] = matrix[r], matrix[i]
        lv = matrix[r][lead]
        matrix[r] = [mrx / float(lv) for mrx in matrix[r]]
        for i in range(rowCount):
            if i != r:
                lv = matrix[i][lead]
                matrix[i] = [iv - lv*rv for rv,iv in zip(matrix[r],matrix[i])]
        lead += 1
    return matrix

def gcd(a, b):
    if a == 0:
        return b
    if b == 0:
        return a
    if a == b:
        return a
    if a > b:
        return gcd(a-b, b)
    return gcd(a, b-a)

def is_prime(n):
    if n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False
    for i in range(3, isqrt(n)+1, 2):
        if n % i == 0:
            return False
    return True

def smoothness_bound(n):
    return int(1.5 * pow(n, 1/4))

def legendre_symbol(a, p):
    ls = pow(a, (p-1)//2, p)
    return ls if ls == 1 else -1

def sieve(n, b):
    primes = [i for i in range(2, b+1) if is_prime(i)]
    factor_base = []
    for p in primes:
        if legendre_symbol(n, p) == 1:
            e = 0
            while n % p == 0:
                n //= p
                e += 1
            factor_base.append((p, e))
    if n > 1:
        return None
    return factor_base

def quad_residue(a, n, p):
    x = pow(a, (p+1)//4, p)
    return x if pow(x, 2, p) == a else n-x

def smooth_sieve(n, b, factor_base):
    sieve = []
    for i in range(1, b+1):
        x = quad_residue(i, n, factor_base[0])
        factors = []
        for p in factor_base:
            e = 0
            while x % p == 0:
                x //= p
                e += 1
            factors.append(e)
        sieve.append((i, factors))
    return sieve


if __name__ == "__main__":
    # test quadratic sieve
    n = 131 * 137
    factors = quadratic_sieve(n)
    print(factors)

    # test quadratic sieve
    n = 131 * 137 * 139
    factors = quadratic_sieve(n)
    print(factors)

    # test quadratic sieve
    n = 131 * 137 * 139 * 149
    factors = quadratic_sieve(n)
    print(factors)

    # test quadratic sieve
    n = 131 * 137 * 139 * 149 * 151
    factors = quadratic_sieve(n)
    print(factors)

    # test quadratic sieve
    n = 131 * 137 * 139 * 149 * 151 * 157
    factors = quadratic_sieve(n)
    print(factors)

    # test quadratic sieve
    n = 131 * 137 * 139 * 149 * 151 * 157 * 163
    factors = quadratic_sieve(n)
    print(factors)

    # test quadratic sieve
    n = 131 * 137 * 139 * 149 * 151 * 157 * 163 * 167
    factors = quadratic_sieve(n)
    print(factors)

    # test quadratic sieve
    n = 131 * 137 * 139 * 149 * 151 * 157 * 163 * 167 * 173
    factors = quadratic_sieve(n)
    print(factors)

    # test quadratic sieve
    n = 131 * 137 * 139 * 149 * 151 * 157 * 163 * 167 * 173 * 179
    factors = quadratic_sieve(n)
    print(factors)

    # test quadratic sieve
    n = 131 * 137 * 139 * 149 * 151 * 157 * 163 * 167 * 173 * 179 * 181
    factors = quadratic_sieve(n)
    print(factors)

    # test quadratic sieve
    n = 131 * 137 * 139 * 149 * 151 * 157 * 163 * 167 * 173 * 179 * 181 * 191
    factors = quadratic_sieve(n)
    print(factors)