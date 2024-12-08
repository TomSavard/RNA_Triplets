import numpy as np
# print(np.__version__)
# We will focus on dissecting the problem into tiny tasks that we can approach and verify easily.
a_sequence = "GGGAAAUUCCUUAAAGGCUUCCGGAU"
theta = 4 # common default value. Theta is the min distance between two bases for them to be considered for pairing.
pair_energy = -1 #to be defined


###########         Utilities


def print_matrix(m, sequence):
    n = len(m)
    # Print the header row with the sequence positions
    print("     ", "   ".join(f"{sequence[i]:>3}" for i in range(n)))  # Adjusted spacing between columns
    for i in range(n):
        # Print each row label with extra space for better alignment
        print(f"{sequence[i]:>3}  ", end="")  # Add more space after the label
        for j in range(n):
            if m[i, j] == 0:
                print("   0  ", end="")  # Extra space for zero values
            else:
                print(f"{m[i, j]:>5.1f}", end=" ")  # Increased space for matrix values
        print()  # Move to the next line for the next row



##### Main functions

def FillMatrix(w) :
    """
    Fills the minimum energy matrix

    Args:
        w (string): RNA sequence

    Returns:
        matrix of size n*n: minimum energy matrix of the sequence
    """
    n = len(w)
    m = np.zeros((n,n))
    for i in range(n): 
        for j in range (i, min(i + theta, n)) :
            m[i,j] = 0
            # this step is unnecessary as we are replacing zeros with zeros but it helps in understanding the algo
    for i in range (n-1,-1,-1): # We start with larger subsequences and iteratively calculate smaller ones.
        for j in range (i + theta + 1, n) :
            # Case A : pos i without partner
            m[i,j] = m[i+1,j]
            # Case B : Pos i and j form a base pair
            if can_pair(w[i], w[j]) :  # if pos i and j can pair 
                m[i,j] = min(m[i,j], m[i+1,j-1] + pair_energy)
            # Case C : pairing with an other base at index k
            for k in range (i + theta + 1, j):
                m[i,j] = min(m[i,j], m[i+1,k-1] + m[k+1,j] + pair_energy)
    return m

def backtracking(i,j,m,w):
    """
    Finds the minimum energy structure

    Args:
        i (int): Begining of the region under consideration
        j (int): End of the region under consideration
        m (matrix n*n): minimum energy matrix
        w (string): RNA sequence of length n

    Returns:
        Something: structure minimizing free energy
    """

    return 0


##### Main

possible_pairs = [("A", "U"), ("G", "C"), ("G", "U"), 
                    ("U", "A"), ("C", "G"), ("U", "G")]

def can_pair(base1,base2):
    return (base1,base2) in possible_pairs

m = FillMatrix(a_sequence)

print("Final matrix:")





# Print the matrix in a formatted manner
print_matrix(m, a_sequence)

