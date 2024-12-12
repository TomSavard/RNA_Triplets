import numpy as np
# print(np.__version__)
# We will focus on dissecting the problem into tiny tasks that we can approach and verify easily.
a_sequence = "GGGAAAUUCCUUA"
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


def can_pair(base1,base2):
    return (base1,base2) in possible_pairs

def print_pairing(w,S):
    for (i,j) in S :
        print(w[i],w[j])

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
    print("size of the sequence: ",n)
    m = np.zeros((n,n))
    for i in range(n): 
        for j in range (i, min(i + theta, n)) :
            m[i,j] = 0
            # this step is unnecessary as we are replacing zeros with zeros but it helps in understanding the algo
    for i in range (n-1,-1,-1): # We start with larger subsequences and iteratively calculate smaller ones.
        print("Processing the index ", i, " corresponding to the letter : ", w[i])
        for j in range (i + theta + 1, n) :
            print("     Based on theta we can consider the index: ", j, " letter: ", w[j])
            # Case A : pos i without partner
            m[i,j] = m[i+1,j]
            # Case B : Pos i and j form a base pair
            if can_pair(w[i], w[j]) :  # if pos i and j can pair
                print( "        we consider the subsequence [", i, ",", j, "]")
                m[i,j] = min(m[i,j], m[i+1,j-1] + pair_energy)
                print("         m[",i,",",j,"] = ",m[i,j] )
            # Case C : pairing with an other base at index k
            for k in range (i + theta + 1, j):
                print("k = ",k)
                m[i,j] = min(m[i,j], m[i+1,k-1] + m[k+1,j] + pair_energy)
        print(m)
    return m

def Backtrack(i,j,m,w,S):
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
    if j-i <= theta :
        return
    else :
        # Case A : Position i left without partner
        if m[i,j] == m[i+1,j] :
            Backtrack(i + 1, j, m, w, S)
            return 
        # Case B : Positions i and j form a base pair 
        if m[i,j] == m[i+1,j-1] + pair_energy :
            S.append((i, j))
            Backtrack(i + 1, j - 1, m, w, S)
            return



            


##### Main

possible_pairs = [("A", "U"), ("G", "C"), ("G", "U"), 
                    ("U", "A"), ("C", "G"), ("U", "G")]

n = len(a_sequence)
m = FillMatrix(a_sequence)
# print("The RNA sequence analysed: ", a_sequence)
# print("Size of the RNA sequence: ", n)

print("minimum energy matrix:")
print_matrix(m, a_sequence)

# S = []
# Backtrack(0,n-1,m,a_sequence,S)
# print("The optimal pairing :", S)

# print_pairing(a_sequence,S)