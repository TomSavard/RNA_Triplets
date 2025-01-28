from global_variables import *


def print_matrix(m, sequence, cell_width=3):
    """
    Prints the matrix with aligned columns by setting a fixed width for each cell.

    Args:
        m (numpy.ndarray): The matrix to print.
        sequence (str): The RNA sequence.
        cell_width (int): The width of each cell (default is 6).
    """
    n = len(m)

    print(" " * cell_width, " ".join(f"{sequence[i]:>{cell_width}}" for i in range(n)))  # Adjust header alignment

    for i in range(n):
        # Print each row label with extra space for better alignment
        print(f"{sequence[i]:>{cell_width}}", end=" ")
        for j in range(n):
            print(f"{int(m[i, j]):>{cell_width}}", end=" ")  # Align matrix values and format as floats with 1 decimal place
        print()  # Move to the next line for the next row





def can_pair(base1,base2):
    return (base1,base2) in possible_pairs





def print_pairing(sequence,S):
    for (i,j) in S :
        print(sequence[i],sequence[j])





def displaySS(pairing, size):
    """
    Display the secondary structure of the RNA sequence

    Args:
        pairing (list): list of pairs of indices of the RNA sequence
        size (int): size of the RNA sequence

    Returns:
        string: secondary structure of the RNA sequence as a parnthesized string

    """
    structure = ["." for _ in range(size)]
    for (i,j) in pairing :
        structure[i] = "("
        structure[j] = ")"
    return "".join(structure)





def parseSS(struct): 
    """
    Parse the secondary structure of the RNA sequence

    Args:
        struct (string): secondary structure of the RNA sequence

    Returns:
        tuple: list of base-pairs and length of the string
    """
    l = []
    n = len(struct)
    last_index = n-1
    for i in range(n):
        if struct[i] == "(":
            for j in range(last_index,i,-1):
                if struct[j] == ")":
                    l.append((i,j))
                    last_index = j-1
                    break
    return (l,n)




if __name__ == "__main__" :
    print("\n--------------------------------------------------")
    print("Test of the utilities\n")

    print("Test of displaySS")
    print(displaySS([(2,8),(3,7)],10)) # Expected output: "..((...))."

    print("\n")

    print("Tets of parseSS")
    print(parseSS("((.(...)))")) # Expected output: ([(0, 9), (1, 8), (3, 7)], 10)
    print("\nEnd of utilities' test")
    print("--------------------------------------------------\n")