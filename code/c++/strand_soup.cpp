#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <ctime>

#include "nussinov.hpp"


// Here we want to implement the strand soup interaction
//This is described in the paper AIMoB_submitted version figure 7

// There are two different parts

// 1. The genneral case (represented with a circle on the figure) with two strands.
// 2. The special case where we have to check for emptiness (represented with a square on the figure).



// The set up 

// We have a set of strands S = {s1, s2, ..., sn} where each strand is a sequence of nucleotides.
// We select a starting strand s and an ending strand r.
// then we have to iteratively compute each element of 2 6D matrixes.
// We are have to go bottom up in the matrices. We follow the lexicographic order (m,s-i,r,j,c)
// Once this will be done, we will have to backtrack to find the optimal secondary structure. (How is the backtrack done ?)




// Fonction pour générer une séquence de nucléotides aléatoire de longueur donnée
std::string generate_random_sequence(int length) {
    const std::string nucleotides = "ACGU";
    std::string sequence;
    for (int i = 0; i < length; ++i) {
        sequence += nucleotides[rand() % 4];
    }
    return sequence;
}



float inf_energy = 1000000; // We should discuss this value. Simply having a positive energy could be sufficient as any other would be negative for the moment. (Temporary solution)
float pair_energy = -1; //! different from the global variable as we need a float type for the std::min
using Matrix6D = std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<float>>>>>>;

Matrix6D initialize_matrix(int s_size, int r_size, int i_size, int j_size, int m_size, int c_size) {
    return Matrix6D(s_size, std::vector<std::vector<std::vector<std::vector<std::vector<float>>>>>(
        r_size, std::vector<std::vector<std::vector<std::vector<float>>>>(
        i_size, std::vector<std::vector<std::vector<float>>>(
        j_size, std::vector<std::vector<float>>(
        m_size, std::vector<float>(
        c_size, 0))))));
}


//! I need to add a list of boolean to keep track of which strand is still available. 



//// We will use a recursive implementation. 
//! The implementation is iterative not recursive the dynamic scheme lies in how the iteration is handled and the interdependance of each iteration.
//! For the bottom up implementation. We use the lexicographic order described in the paper. (m,s,-i,r,j,c) where m is the number of strands left in the soup, s is the index of the first strand, i is the starting index in the first strand, r is the index of the second strand, j is the ending index in the second strand and c is the connectivity bit.

void AuxiliaryMatrix(int starting_strand, int final_strand, std::unordered_map<int, std::string> strands, Matrix6D& M_bar){
    /**
     * @brief 
     * 
     * Detailled description : TO BE DONE
     * 
     * @param starting_strand : index of the first strand inside strands
     * @param final_strand : index of the second strand inside strands
     * @param strands : the dictionary of strands (int = index, string = sequence)
     * @param M_bar : the auxiliary matrix to fill
     * 
     * @return void
     */

    // ! disjunction based on the emptiness of the strands
    int m_max = strands.size();
    if (m_max <= 2){m_max = 0;}else {m_max = m_max - 2;}
    int s_max = strands.size();
    int i_max = strands[starting_strand].length(); // Here we suppose each strand of same length but we should modify the loop each time
    int r_max = strands.size();    
    int j_max = strands[final_strand].length(); // Here we suppose each strand of same length but we should modify the loop each time

    for (int m=0; m < m_max; m++){
        for (int s=0; s < s_max; s++){
            for (int i=strands[s].length(); i >= 0; i--){ // We need to start from the end of the strand
                for (int r=0; r < r_max; r++){
                    for (int j=0; j < strands[r].length(); j++){
                        for (int c=0; c<2; c++){
                            // Computation of M[m,s,i,r,j,c]
                            if (i > strands[s].length()) { // if s empty
                                if (c == 1){// disconnected strands : energy = +inf
                                    M_bar[m][s][i][r][j][c] = inf_energy;
                                }
                                else{
                                    if (m == 0){// we call nussinov on r[1,j]
                                        M_bar[m][s][i][r][j][c] = Nussinov(strands[r],1,j);
                                    }
                                    else {// We insert a strand t before the strand r and we minimize over all possible t
                                        float temp;
                                        for (int t = 0; t < m-1; t++){
                                            // Todo : We need to check if the strand t is available. 
                                            temp = std::min(temp,M_bar[m-1][t][1][r][j][1]);
                                        }
                                        M_bar[m][s][i][r][j][c] = temp;
                                    }
                                }
                        
                            }
                            else { // if s not empty
                                if (j> strands[r].length()){// if r empty
                                    if (c == 1){// disconnected strands : energy = +inf
                                        M_bar[m][s][i][r][j][c] = inf_energy; //Todo : inf to be modified
                                    }
                                    else{
                                        if (m == 0){// we call nussinov on s
                                            M_bar[m][s][i][r][j][c] = Nussinov(strands[s],i,strands[s].length()-1); //! careful on the final index. It is written |s| on the paper so it should be |s| -1 for us ?
                                        }
                                        else {// We insert a strand t after the strand s and we minimize over all possible t
                                            float temp;
                                            for (int t = 0; t < m-1; t++){
                                                // Todo : We need to check if the strand t is available. 
                                                temp = std::min(temp,M_bar[m-1][s][i][t][strands[t].length()-1][1]);
                                            }
                                            M_bar[m][s][i][r][j][c] = temp;
                                        }
                                    }
                        
                                }
                                else {
                                    // Todo : We need the value of the main matrix M[m][s][i][r][j][c] but it isn't computed yet :/
                                    // ! Do we need to compute both matrices in parallel ? 
                                }
                            }

                        }
                    }
                }
            }
        }
    }

}


void MainMatrix(int starting_strand, int final_strand, std::unordered_map<int, std::string> strands, Matrix6D& M){
    /**
     * @brief 
     * 
     * Detailled description : TO BE DONE
     * 
     * @param starting_strand : index of the first strand inside strands
     * @param final_strand : index of the second strand inside strands
     * @param strands : the dictionary of strands (int = index, string = sequence)
     * @param M : the main matrix to fill
     * 
     * @return void
     */

    int m_max = strands.size();
    if (m_max <= 2){m_max = 0;}else {m_max = m_max - 2;}
    int s_max = strands.size();
    int i_max = strands[starting_strand].length(); // Here we suppose each strand of same length but we should modify the loop each time
    int r_max = strands.size();    
    int j_max = strands[final_strand].length(); // Here we suppose each strand of same length but we should modify the loop each time

    // Iteration over the matrix in lexicographic order
    for (int m=0; m < m_max; m++){
        for (int s=0; s < s_max; s++){
            for (int i=strands[s].length(); i >= 0; i--){ // We need to start from the end of the strand
                for (int r=0; r < r_max; r++){
                    for (int j=0; j < strands[r].length(); j++){
                        for (int c=0; c<2; c++){
                            // 4 cases to consider and take the min

                            // Case 1 : i left unpaired. We need to check if the strand is empty
                            float ca@se1 = M_bar[m][s][i+1][r][j][c]; // Todo : here we might exceed the index of s in the matrix. We need to find a solution for that.
                            // Case 2 : i is paired with k belonging to s. We call Nussinov on (s,i+1,k-1) and need M_bar[m][s][k+1][r][j][c]
                            // We have to minimize over k, 
                            float case2 = inf_energy;
                            for (int k = i+1; k < strands[s].length(); k++){
                                case2 = std::min(case2, Nussinov(strands[s],i+1,k-1) + M_bar[m][s][k+1][r][j][c] + pair_energy); //Todo : pair_energy = E(s,i,r,j) we can improve the code and be more specific.
                            }
                            // Case 3 : i is paired with k belonging to t. We need M_bar[m'][s][i+1][t][k-1][0] AND M_bar[m-m'-1][t][k+1][r][j][c]
                            // We have to minimize over t, k and m
                            float case3;
                            for (int t = 0; t< m; t++){
                                for (int k = 0; k< strands[t].length(); k++){
                                    for (int m_prime = 0; m_prime <= m-1 ; m_prime++){
                                        case3 = std::min(case3, pair_energy + M_bar[m_prime][s][i+1][t][k-1][0] + M_bar[m-m_prime-1][t][k+1][r][j][c]);
                                    }
                                }
                            }
                            // Case 4 : i is paired with k belonging to r. We call Nussinov on (r,k+1,j) and need M_bar[m][s][i+1][r][k-1][0]
                            float case_4;
                            for (int k = 0; k < strands[r].length(); k++){
                                case_4 = std::min(case_4, pair_energy + M_bar[m][s][i+1][r][k-1][c] + Nussinov(strands[r],k+1,j));
                            }
                            M[m][s][i][r][j][c] = std::min({case1,case2,case3,case_4});
                        }
                    }
                }
            }
        }
    }
    
}


int main() {
    srand(time(0)); // Initialiser le générateur de nombres aléatoires

    // Générer un ensemble de séquences
    int m = 2; // Nombre de séquences à générer
    int sequence_length = 10; // Longueur de chaque séquence
    std::vector<std::string> sequences;

    for (int i = 0; i < m; ++i) {
        sequences.push_back(generate_random_sequence(sequence_length));
    }

    // Afficher les éléments
    std::cout << "  Number of strands : m = " << m << std::endl;
    std::cout << "  Generated sequences:" << std::endl;
    for (const auto& seq : sequences) {
        std::cout <<"   " <<  seq << std::endl;
    }
    std::string a_sequence = sequences[0];
    std::vector<std::vector<int>> energy_matrix;
    FillMatrix(a_sequence, energy_matrix);
    print_matrix(energy_matrix,a_sequence);




    // Test accès matrice 6D

    int s_size = 3;
    int r_size = 3;
    int i_size = 10;
    int j_size = 10;
    int m_size = 3;
    int c_size = 2; // Connectivité : 0 ou 1
    Matrix6D M = initialize_matrix(s_size, r_size, i_size, j_size, m_size, c_size);

    // Afficher les éléments du vecteur de la dernière dimension
    std::cout << "M[" << 0 << "][" << 0 << "][" << 0 << "][" << 0 << "][" << 0 << "] = ";
    for (int k = 0; k < c_size; ++k) {
        std::cout << M[0][0][0][0][0][k] << " ";
    }
    std::cout << std::endl;


    // Test du dictionnaire
    std::unordered_map<int, std::string> strands;
    for (int i =0; i < m ; i++){
        strands[i] = generate_random_sequence(sequence_length);
    }
    std::cout << "Generated strands:" << std::endl;
    for (const auto& pair : strands){
        std::cout << "Index " << pair.first << " : " << pair.second << std::endl;
    }



    return 0;
}