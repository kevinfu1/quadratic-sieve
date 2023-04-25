from math import isqrt
from itertools import chain

def gcd(a, b):
    """
    This function computes the greatest common divisor (GCD) of two numbers using the Euclidean algorithm.

    :param a: The first number.
    :param b: The second number.
    :return: The GCD of a and b.

    Time complexity: O(log(min(a, b)))
    """

    if b == 0:  # If b is 0, the GCD is a
        return a
    elif a >= b:  # If a is greater than or equal to b, continue with the Euclidean algorithm
        return gcd(b, a % b)  # Recursively call the function with b and a mod b
    else:  # If b is greater than a, swap the values and call the function again
        return gcd(b, a)

def transpose(matrix: list) -> list:
    """
    Transposes a matrix.

    :param matrix: A list of lists representing the matrix to be transposed.
    :return: A new matrix with the rows and columns of the original matrix swapped.

    Time Complexity: O(r*c), where r and c are the number of rows and columns respectively
    """

    new_matrix = []
    for i in range(len(matrix[0])):  # Iterate over each column of the matrix. O(c)
        new_row = []
        for row in matrix:  # Iterate over each row of the matrix. O(r)
            new_row.append(row[i])  # Append the ith element of each row to a new row. O(1)
        new_matrix.append(new_row)  # Append the new row to the new matrix. O(1)
    return new_matrix  # O(1)


#reduced form of gaussian elimination, finds rref and reads off the nullspace
#https://www.cs.umd.edu/~gasarch/TOPICS/factoring/fastgauss.pdf

def gauss_elim(M):

    #M = optimize(M)
    marks = [False]*len(M[0])
    
    for i in range(len(M)): #do for all rows
        row = M[i]
        #print(row)
        
        for num in row: #search for pivot
            if num == 1:
                #print("found pivot at column " + str(row.index(num)+1))
                j = row.index(num) # column index
                marks[j] = True
                
                for k in chain(range(0,i),range(i+1,len(M))): #search for other 1s in the same column
                    if M[k][j] == 1:
                        for i in range(len(M[k])):
                            M[k][i] = (M[k][i] + row[i])%2
                break
    
    M = transpose(M)
        
    sol_rows = []
    for i in range(len(marks)): #find free columns (which have now become rows)
        if marks[i]== False:
            free_row = [M[i],i]
            sol_rows.append(free_row)
    
    if not sol_rows:
        return("No solution found. Need more smooth numbers.")
    print("Found {} potential solutions".format(len(sol_rows)))
    return sol_rows,marks, M

def solve_row(sol_rows,M,marks,K=0):
    solution_vec, indices = [],[]
    free_row = sol_rows[K][0] # may be multiple K
    for i in range(len(free_row)):
        if free_row[i] == 1: 
            indices.append(i)
    for r in range(len(M)): #rows with 1 in the same column will be dependent
        for i in indices:
            if M[r][i] == 1 and marks[r]:
                solution_vec.append(r)
                break
            
    solution_vec.append(sol_rows[K][1])       
    return(solution_vec)
    
def solve(solution_vec,smooth_nums,xlist,N):
    
    solution_nums = [smooth_nums[i] for i in solution_vec]
    x_nums = [xlist[i] for i in solution_vec]
    
    Asquare = 1
    for n in solution_nums:
        Asquare *= n
        
    b = 1
    for n in x_nums:
        b *= n

    a = isqrt(Asquare)
    
    factor = gcd(b-a,N)
    return factor
