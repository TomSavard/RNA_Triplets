#! Organization :

# 1. Dependencies : contains the libraries used in the code
# 2. Main functions : contains the main functions of the code
# 3. Main : code where we execute the main functions and initialize the variables


#! Il faut absolument gérer les paramètres de façons rigoureuse et se concentrer sur les test unitaires !!!!
# ! J'ai passé deux jours à vérifier mon implémentation du backtracking et regarder des videos pour au final me rendre compte que mon paramètre theta n'était pas le même entre le FillMatrix et le BackTrack...

###########         Dependencies

import numpy as np
# print(np.__version__)
import subprocess

from utilities import *
from global_variables import *

import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from PIL import Image







###########         Main functions



def FillMatrix(sequence, theta = 3, pair_energy = -1):
    """
    Fills the minimum energy matrix

    Args:
        sequence (string): RNA sequence

    Returns:
        matrix of size n*n: minimum energy matrix of the sequence
    """
    print("\n-------------------")
    print("Start of FillMatrix")

    len_seq = len(sequence)
    # ----------- Initialisation of the matrix with zeros -----------
    m = np.zeros((len_seq,len_seq))
    for i in range (len_seq): 
        for j in range (i, min(i + theta, len_seq)) :
            m[i,j] = 0 # this step is unnecessary as we are replacing zeros with zeros but it helps in understanding the algo

    # ----------- Filling the matrix -----------
    for i in range (len_seq-1,-1,-1): # We start with larger subsequences and iteratively calculate smaller ones.
        for j in range (i + theta + 1, len_seq) :
            # Case A : pos i without partner
            m[i,j] = m[i+1,j]

            # Case B : Pos i and j form a base pair
            if can_pair(sequence[i], sequence[j]) :
                m[i,j] = min(m[i,j], m[i+1,j-1] + pair_energy)

            # Case C : pairing with an other base at index k
            for k in range (i + theta + 1, j):
                if can_pair(sequence[i], sequence[k]) :
                    # print("k = ",k)
                    m[i,j] = min(m[i,j], m[i+1,k-1] + m[k+1,j] + pair_energy)

    print("End of FillMatrix")
    print("-------------------\n")
    return m







def Backtrack(i, j, m, sequence, S, theta = 3):
    """
    Finds the minimum energy structure

    Args:
        i (int): Begining of the region under consideration
        j (int): End of the region under consideration
        m (matrix n*n): minimum energy matrix
        sequence (string): RNA sequence of length n

    Returns:
        Something: structure minimizing free energy
    """

    # print("\n-------------------")
    # print("Start of Backtrack")
    # print("i = ",i, " j = ",j)
    # print(sequence[i] + "       " + sequence[j])

    if j-i <= theta :
        # print("CaseD")
        return
    else :
        # Case A : Position i unpaired
        if m[i,j] == m[i+1,j] :
            # print("CaseA")
            Backtrack(i + 1, j, m, sequence, S)
            return 
        # elif m[i,j] == m[i,j-1]: # It also works but gives an other possible structure
        #     Backtrack(i, j-1, m, sequence, S)
        # Case B : Positions i and j form a base pair 
        elif m[i,j] == m[i+1,j-1] + pair_energy :
            # print("CaseB")
            S.append((i, j))
            Backtrack(i + 1, j - 1, m, sequence, S)
            return
        # Case C : Position i pairs with an other base k < j
        else :
            for k in range(i + theta + 1, j) :
                if  m[i,j] == m[i+1,k-1] + m[k+1,j] + pair_energy :
                    print("CaseC")
                    S.append((i, k))
                    Backtrack(i+1, k-1, m, sequence, S)
                    Backtrack(k + 1, j, m, sequence, S)
                    return
        
    




            


###########         Main

if __name__ == "__main__":

    a_sequence = "GCAACUGGCAC"
    len_seq = len(a_sequence)

    print("The RNA sequence : ", a_sequence)
    print("Size of the RNA sequence: ", len_seq)

    m = FillMatrix(a_sequence)


    print("Minimum energy matrix:")
    print_matrix(m, a_sequence)

    S = []
    Backtrack(0, len_seq-1, m, a_sequence, S, theta = 3)
    print("The optimal pairing :", S)
    structure = displaySS(S,len_seq)
    print(structure)


    rna_input = f"{a_sequence}\n{structure}\n"
    with open("../../data/rna_input.txt", "w") as f:
        f.write(rna_input)

    subprocess.run(["RNAplot", "--output-format=svg", "../../data/rna_input.txt"])
    print("Visualization saved as rna.svg")


    # svg_file = "rna.svg"
    # output_png = "rna_input_ss.png"
    # Image.open(svg_file).save(output_png)

    # # Display the PNG version of the SVG
    # img = plt.imread(output_png)
    # fig, ax = plt.subplots(figsize=(8, 8))
    # ax.imshow(img)
    # ax.axis("off")  # Hide axes
    # plt.show()

    from IPython.display import SVG, display

    # Display the SVG
    svg_file = "../../data/rna.svg"
    display(SVG(svg_file))