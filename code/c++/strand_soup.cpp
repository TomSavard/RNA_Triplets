#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <ctime>

#include "nussinov.hpp"


//=============================================================================================//


// Here we want to implement the strand soup interaction
//This is described in the paper AIMoB_submitted version figure 7

// There are two different parts

// 1. The genneral case (represented with a circle on the figure) with two strands.
// 2. The special case where we have to check for emptiness (represented with a square on the figure).



// The set up 

// We have a set of strands S = {s1, s2, ..., sn} where each strand is a sequence of nucleotides.
// We select a starting strand s and an ending strand r.
// then we have to iteratively compute each element of a 6D matrixes.
// We are have to go bottom up in the matrices. We follow the lexicographic order (m,s-i,r,j,c)
// Once this will be done, we will have to backtrack to find the optimal secondary structure. (How is the backtrack done ?)



//=============================================================================================//


// Function to generate a random sequence of nucleotides of a given length 
std::string generate_random_sequence(int length) {
    /**
     * @brief generate a random sequence of nucleotides of a given length
     * 
     * @param length (int): the length of the sequence
     * 
     * @return std::string : the generated sequence
     */
    const std::string nucleotides = "ACGU";
    std::string sequence;
    for (int i = 0; i < length; ++i) {
        sequence += nucleotides[rand() % 4];
    }
    return sequence;
}


class Matrix6D{
    public:
    Matrix6D(int m_size, int s_size, int i_size, int r_size, int j_size, int c_size)
        : m_size(m_size), s_size(s_size), i_size(i_size), r_size(r_size), j_size(j_size), c_size(c_size) {
        data = std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<float>>>>>>(
            m_size, std::vector<std::vector<std::vector<std::vector<std::vector<float>>>>>(
            s_size, std::vector<std::vector<std::vector<std::vector<float>>>>(
            i_size, std::vector<std::vector<std::vector<float>>>(
            r_size, std::vector<std::vector<float>>(
            j_size, std::vector<float>(
            c_size, 0))))));
    }
    float& operator()(int m, int s, int i, int r, int j, int c) {
        return data[m][s][i][r][j][c];
    }

    int get_m_size() const { return m_size; }
    int get_s_size() const { return s_size; }
    int get_i_size() const { return i_size; }
    int get_r_size() const { return r_size; }
    int get_j_size() const { return j_size; }
    int get_c_size() const { return c_size; }

private:
    int m_size, s_size, i_size, r_size, j_size, c_size;
    std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<float>>>>>> data;
};




float inf_energy = 1000000; // We should discuss this value. Simply having a positive energy could be sufficient as any other would be negative for the moment. (Temporary solution)




float GeneralCaseMinimization (int m, int s, int i, int r, int j, int c, std::unordered_map<int, std::string> strands, Matrix6D& M,  std::unordered_map<std::string, std::vector<std::vector<float>>> Nussinov_matrices){
    /**
     * @brief Computes the minimum energy for the general case.
     * 
     * Detailled description : TO BE DONE
     * 
     * @param m : the number of remaining strands in the soup
     * @param s : the index of the starting strand
     * @param i : the index of the starting nucleotide of the first strand
     * @param r : the index of the ending strand
     * @param j : the index of the ending nucleotide of the last second strand
     * @param c : the connectivity bit
     * @param strands : the dictionary of strands (int = index, string = sequence)
     * @param M : the energy_matrix
     * 
     * @return float : the minimum energy
     */

    float min_value = inf_energy;

    // Case 1 : i is left unpaired //! we have top deal with the border effect. We must add a column to the indices i and j. This column at index strands.length() will be the empty case. // We can also have a border effect in front (for the second strand)
    min_value = M(m,s,i+1,r,j,c);

    // Case 2 : i is paired to k of strand s
    for (int k=i+1; k<int(strands.at(s).length()); k++){
        if (can_pair(strands.at(s)[i],strands.at(s)[k])){
            float value = pair_energy + Nussinov_matrices[strands.at(s)][i+1][k-1] + M(m,s,k+1,r,j,c);
            if (value < min_value){
                min_value = value;
            }
        }

    }

    // Case 3 : i is paired to k of new strand t
    if (m>1){
        for (int t=0; t<M.get_s_size(); t++){
            for (int m1=0; m1<m; m1++){
                int m2 = m - m1 - 1;
                for (int k=0; k<int(strands.at(t).length()); k++){
                    if (can_pair(strands.at(s)[i],strands.at(t)[k])){
                        float value = pair_energy + M(m1,s,i+1,t,k-1,0) + M(m2,t,k+1,r,j,c);
                        if (value < min_value){
                            min_value = value;
                        }
                    }
                }
            }
        }   
    }

    // Case 4 : i is paired to k of strand r

    for (int k=0; k<int(strands.at(r).length()); k++){
        if (can_pair(strands.at(s)[i],strands.at(r)[k])){
            float value = pair_energy + M(m,s,i+1,r,k-1,0) + Nussinov_matrices[strands.at(r)][k+1][j];
            if (value < min_value){
                min_value = value;
            }
        }
    }

    return min_value;

}







void MainAuxiliaryMatrix( std::unordered_map<int, std::string> strands, Matrix6D& M, std::unordered_map<std::string, std::vector<std::vector<float>>> Nussinov_matrices){
    /**
     * @brief Computes the 6D matrix of the strand soup.
     * 
     * Detailled description : Sebastian Will suggested the idea of only computing the "square" 6D matrix and replacing the first one by a function with a disjunction of cases.
     * 
     * int starting_strand, int final_strand // no longer using these arguments
     * @param starting_strand : index of the first strand inside strands
     * @param final_strand : index of the second strand inside strands
     * @param strands : the dictionary of strands (int = index, string = sequence)
     * @param M : the main matrix to fill
     * 
     * @return void
     */

    std::cout << "MainAuxiliaryMatrix" << std::endl;
    // ================== Setting the range for the iterations ================== //
    int m_max = M.get_m_size();
    if (m_max <= 2){m_max = 0;}else {m_max = m_max - 2;}
    int s_max = M.get_s_size();
    int i_max = M.get_i_size(); // Here we suppose each strand of same length but we should modify the loop each time
    int r_max = M.get_r_size();
    int j_max = M.get_j_size(); // Here we suppose each strand of same length but we should modify the loop each time
    int c_max = M.get_c_size();


    // ================== Lexicographic Order (m,s,-i,r,j,c) ================== //
    for (int m=0; m<=m_max; m++){
        std::cout << "m = " << m << std::endl;
        for (int s=0; s< s_max; s++){
            for (int i=i_max-1; i>=0; i--){
                for (int r=0; r<r_max; r++){
                    for (int j=0; j<j_max; j++){
                        for (int c=0; c<c_max; c++){
                            std::cout << "M[" << m << "][" << s << "][" << i << "][" << r << "][" << j << "][" << c << "] = " << M(m,s,i,r,j,c) << std::endl;
                            if (i >= int(strands[s].length())){// if s is empty //! How can it be empty since i < i_max ?
                                if(c==1){M(m,s,i,r,j,c) = inf_energy;}
                                else{
                                    if (m==0){M(m,s,i,r,j,c) = Nussinov_matrices[strands.at(r)][0][j];}
                                    else{// we select the min M[m-1][t][0][r][j-1][1] over all t
                                        float min_value = inf_energy;
                                        for (int t=0; t<s_max; t++){
                                            if (M(m-1,t,0,r,j,c) < min_value){min_value = M(m-1,t,0,r,j,c);}}
                                        M(m,s,i,r,j,c) = min_value;
                                    }
                                }
                            }
                            else{
                                if (j>=int(strands[r].length())){//if r is empty
                                    if (c==1){M(m,s,i,r,j,c) = inf_energy;}
                                    else {
                                        if (m==0){M(m,s,i,r,j,c) = Nussinov_matrices[strands.at(s)][i][strands.at(s).length()];}
                                        else{// we select the min M[m-1][s][i][t][size(t)-1][1] over all t
                                            float min_value = inf_energy;
                                            for (int t=0; t<s_max; t++){
                                                if (M(m-1,s,i,t,strands[t].length()-1,c) < min_value){min_value = M(m-1,s,i,t,strands[t].length()-1,c);}}
                                            M(m,s,i,r,j,c) = min_value;
                                        }
                                    }
                                }
                                else {
                                    GeneralCaseMinimization(m,s,i,r,j,c,strands,M,Nussinov_matrices);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}


int main() {
    std::cout << "\n==================== Strand Soup ====================" << std::endl;
    srand(time(0)); // Initialiser le générateur de nombres aléatoires

    //============= Setting the parameters =============//
    std::cout << "\n========== Setting the parameters ==========" << std::endl;
    int m = 2; // Number of sequences to generate
    int sequence_length = 5; // length of the sequences
    std::cout << "  Number of strands : m = " << m << std::endl;
    std::cout << "  Length of the sequences : " << sequence_length << std::endl;


    //============= Generation of the strands =============//
    std::cout << "\n========== Generation of the strands ==========" << std::endl;
    std::unordered_map<int, std::string> strands;
    for (int i =0; i < m ; i++){
        strands[i] = generate_random_sequence(sequence_length);
    }
    std::cout << "Generated strands:" << std::endl;
    for (const auto& pair : strands){
        std::cout << "  Index " << pair.first << " : " << pair.second << std::endl;
    }


    //============= Nussinov matrices =============//
    std::cout << "\n========== Nussinov matrices ==========" << std::endl;
    std::unordered_map<std::string, std::vector<std::vector<float>>> Nussinov_matrices;
    for (const auto& pair : strands) {
        const auto& seq = pair.second;
        std::cout << "  Energy matrix for sequence : " << seq << std::endl;
        std::vector<std::vector<float>> energy_matrix;
        FillMatrix(seq, energy_matrix);
        Nussinov_matrices[seq] = energy_matrix;
        print_matrix(energy_matrix,seq);
        std::cout << std::endl;
    }



    //============= Test of the 6D matrix =============//
    std::cout << "========== Test of the 6D matrix ==========" << std::endl;
    int m_size = m; // Total number of strands we want to use
    int s_size = m; // Total number of different strands
    int i_size = sequence_length; //Max length of a strand
    int r_size = s_size; // We have the same set of strands as for s
    int j_size = sequence_length; // We have the same length for the strands
    int c_size = 2; // Connectivité : 0 ou 1
    Matrix6D M(m_size, s_size, i_size, r_size, j_size, c_size);

    // Afficher les éléments du vecteur de la dernière dimension
    std::cout << "M[" << 0 << "][" << 0 << "][" << 0 << "][" << 0 << "][" << 0 << "] = ";
    for (int k =0; k < M.get_c_size();k++){
        std::cout << M(0, 0, 0, 0, 0, k) << " ";
    }
    std::cout << std::endl;


    //============= Test of MainAuxiliaryMatrix =============//
    std::cout << "========== Test of MainAuxiliaryMatrix ==========" << std::endl;
    MainAuxiliaryMatrix(strands, M, Nussinov_matrices);



    return 0;
}