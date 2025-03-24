#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <ctime>
#include <limits>

#include "nussinov.hpp"
#include "utilities.hpp"


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





class Matrix6D{
    public:
    Matrix6D(int m_size, int s_size, int i_size, int r_size, int j_size, int c_size)
        : m_size(m_size), s_size(s_size), i_size(i_size), r_size(r_size), j_size(j_size), c_size(c_size) {
        data = std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<float>>>>>>(
            m_size, std::vector<std::vector<std::vector<std::vector<std::vector<float>>>>>(
            s_size, std::vector<std::vector<std::vector<std::vector<float>>>>(
            i_size + 2, std::vector<std::vector<std::vector<float>>>( //! we add one column before and after for border effect
            r_size, std::vector<std::vector<float>>(
            j_size + 2, std::vector<float>( //! we add one column before and after for border effect
            c_size, 0))))));
    }
    // Method to access the elements of the matrix.
    float& operator()(int m, int s, int i, int r, int j, int c) {
        if (s <= 0 || s > s_size) {
            throw std::out_of_range("Index s is out of range");
        }
        if (r <= 0 || r > r_size) {
            throw std::out_of_range("Index r is out of range");
        }
        if (i < 0 || i >= i_size + 2) {
            throw std::out_of_range("Index i is out of range");
        }
        if (j < 0 || j >= j_size + 2) {
            throw std::out_of_range("Index j is out of range");
        }
        return data[m][s-1][i][r-1][j][c]; // 1-based indexing for the strands and the bases. (The 0th column is for the border effect)
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






const float inf_energy = std::numeric_limits<float>::infinity();






float GeneralCaseMinimization (int m, int s, int i, int r, int j, int c, std::unordered_map<int, std::string> strands, Matrix6D& M,  std::unordered_map<std::string, Matrix2D> Nussinov_matrices){
    /**
     * @brief Computes the minimum energy for the general case.
     * 
     * Detailled description : TO BE DONE
     * 
     * @param m : the number of remaining strands in the soup
     * @param s : the index of the starting strand (1-based)
     * @param i : the index of the starting nucleotide of the first strand (1-based
     * @param r : the index of the ending strand (1-based)
     * @param j : the index of the ending nucleotide of the last second strand (1-based)
     * @param c : the connectivity bit
     * @param strands : the dictionary of strands (int = index, string = sequence)
     * @param M : the energy_matrix
     * 
     * @return float : the minimum energy
     */

    //! we have added a column at the end of i and at the begining of j. We have to shift the indices of j by -1 when accessible the bases of the strands (no pb because if j==0, the function isn't called)
    float min_value = inf_energy;

    // Case 1 : i is left unpaired
    min_value = M(m,s,i+1,r,j,c);

    // Case 2 : i is paired to k of strand s
    for (int k=i+1; k<=int(strands.at(s).length())-1; k++){
        if (can_pair(strands.at(s)[i],strands.at(s)[k])){
            float value = pair_energy + Nussinov_matrices[strands.at(s)](i+1,k-1) + M(m,s,k+1,r,j,c);
            if (value < min_value){
                min_value = value;
            }
        }

    }

    // Case 3 : i is paired to k of new strand t
    if (m>=1){
        for (int t=1; t<=M.get_s_size(); t++){
            for (int m1=0; m1<m; m1++){
                int m2 = m - m1 - 1;
                for (int k=1; k<=int(strands.at(t).length())-1; k++){
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

    for (int k=1; k<=int(strands.at(r).length())-1; k++){
        if (can_pair(strands.at(s)[i],strands.at(r)[k])){
            if (k == int(strands.at(r).length()-1)){
                float value = pair_energy + M(m,s,i+1,r,k-1,0);
                if (value < min_value){
                    min_value = value;
                }}
            else {
                float value = pair_energy + M(m,s,i+1,r,k-1,0) + Nussinov_matrices[strands.at(r)](k+1,j);
                if (value < min_value){min_value = value;
                }}
        }
    }
    return min_value;
}






void MainAuxiliaryMatrix( std::unordered_map<int, std::string> strands, Matrix6D& M, std::unordered_map<std::string, Matrix2D> Nussinov_matrices){
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

    // ================== Lexicographic Order (m,s,-i,r,j,c) ================== //
    for (int m=0; m<M.get_m_size(); m++){
        for (int s=1; s<= M.get_s_size(); s++){
            for (int i= M.get_i_size()+1; i>=0; i--){ //! we have added one column before and after for border effect
                for (int r=1; r<=M.get_r_size(); r++){
                    for (int j=0; j<=M.get_j_size()+1; j++){ //! we have added one column before and after for border effect
                        for (int c=0; c<=1; c++){
                            // std::cout << "M(" << m << "," << s << "," << i << "," << r << "," << j << "," << c << ")" << std::endl;
                            if (i > int(strands.at(s).length()-1)){ //! this is the case where the strand s is empty (border effect). Be careful the strands.length must be decreased by one because of the $ sign at the beginning of the strands!
                                // std::cout << "CASE A" << std::endl;
                                if(c==1){
                                    M(m,s,i,r,j,c) = inf_energy;}
                                else{
                                    if (m==0){
                                        if (j == 0) {
                                            M(m,s,i,r,j,c) = 0;} //additional condition in case both are empty
                                        else {
                                            // std::cout << "Nussinov  " << strands.at(r) << " " << 0 << " " << j;
                                            M(m,s,i,r,j,c) = Nussinov_matrices[strands.at(r)](1,j);}
                                    }
                                    else{// we select the min M[m-1][t][0][r][j][1] over all t
                                        float min_value = inf_energy;
                                        for (int t=1; t<=M.get_s_size(); t++){
                                            if (M(m-1,t,1,r,j,1) < min_value){min_value = M(m-1,t,1,r,j,1);}}
                                        M(m,s,i,r,j,c) = min_value;
                                    }
                                }
                            }
                            else{
                                // std::cout << "CASE B" << std::endl;
                                if (j<1){//if r is empty
                                    if (c==1){
                                        M(m,s,i,r,j,c) = inf_energy;}
                                    else {
                                        if (m==0){
                                            if (i ==0){
                                                M(m,s,i,r,j,c) = 0;}
                                            else{
                                                M(m,s,i,r,j,c) = Nussinov_matrices[strands.at(s)](i,strands.at(s).length()-1);}
                                            }
                                        else{// we select the min M[m-1][s][i][t][size(t)-1][1] over all t
                                            float min_value = inf_energy;
                                            for (int t=1; t<=M.get_s_size(); t++){
                                                if (M(m-1,s,i,t,strands[t].length()-1,1) < min_value){min_value = M(m-1,s,i,t,strands[t].length()-1,1);}} //! again be careful of the length of the strands real length = length -1 (-1 because of the $ sign).
                                            M(m,s,i,r,j,c) = min_value;
                                        }
                                    }
                                }
                                else {
                                    // std::cout << "CASE C" << std::endl;
                                    M(m,s,i,r,j,c) = GeneralCaseMinimization(m,s,i,r,j,c,strands,M,Nussinov_matrices);
                                }
                                
                            }
                            std::cout << "M(" << m << "," << s << "," << i << "," << r << "," << j << "," << c << ") = " << M(m,s,i,r,j,c) << std::endl;
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
    int m_start = 2; // Number of sequences to generate
    int sequence_length = 3; // length of the sequences
    std::cout << "  Number of strands : m = " << m_start << std::endl;
    std::cout << "  Length of the sequences : " << sequence_length << std::endl;
    std::cout << "  Theta : " << theta << std::endl;
    std::cout << "  Pair energy : " << pair_energy << std::endl;


    //============= Generation of the strands =============//
    std::cout << "\n========== Generation of the strands ==========" << std::endl;
    std::unordered_map<int, std::string> strands;
    // for (int i =0; i < m ; i++){
    //     strands[i] = generate_random_sequence(sequence_length);
    // }
    strands[1] = "$UUA";
    strands[2] = "$AAC";
    std::cout << "Generated strands:" << std::endl;
    for (const auto& pair : strands){
        std::cout << "  Index " << pair.first << " : " << pair.second << std::endl;
    }



    //============= Nussinov matrices =============//
    std::cout << "\n========== Nussinov matrices ==========" << std::endl;
    std::unordered_map<std::string, Matrix2D> Nussinov_matrices;
    for (const auto& pair : strands) {
        const auto& seq = pair.second;
        std::cout << "  Energy matrix for sequence : " << seq << std::endl;
        Matrix2D energy_matrix(sequence_length, sequence_length);
        FillMatrix(seq, energy_matrix);
        Nussinov_matrices[seq] = energy_matrix;
        print_matrix(energy_matrix,seq);
        std::cout << std::endl;
    }



    //============= Test of the 6D matrix =============//
    std::cout << "========== Test of the 6D matrix ==========" << std::endl;
    int m_size = m_start - 1 ; // Total number of dimensions for m (if we have 2 strands we only have the case m = 0 which is m_size = 1)
    int s_size = m_start; // Total number of different strands
    int i_size = sequence_length; //Max length of a strand
    int r_size = m_start; // We have the same set of strands as for s
    int j_size = sequence_length; // We have the same length for the strands
    int c_size = 2; // Connectivité : 0 ou 1
    Matrix6D M(m_size, s_size, i_size, r_size, j_size, c_size);

    // // print the empty matrix
    // for (int m = 0; m < m_size; m++){
    //     for (int s = 1; s <= s_size; s++){
    //         for (int i = i_size+1; i >=0; i--){
    //             for (int r = 1; r <= r_size; r++){
    //                 for (int j = 0; j <= j_size+1; j++){
    //                     for (int c = 0; c <= 1; c++){
    //                         std::cout << "M(" << m << "," << s << "," << i << "," << r << "," << j << "," << c << ") = ";
    //                         std::cout << M(m, s, i, r, j, c) << std::endl;
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }


    //============= Test of MainAuxiliaryMatrix =============//
    std::cout << "========== Test of MainAuxiliaryMatrix ==========" << std::endl;
    MainAuxiliaryMatrix(strands, M, Nussinov_matrices);

    std::cout << "========== End of the program ==========" << std::endl;

    std::cout <<" Final result" << std::endl;
    std::cout << "M(0,1,1,2,3,1) = " << M(0,1,1,2,3,1) << std::endl;
    std::cout << "M(0,1,1,1,3,1) = " << M(0,1,1,1,3,1) << std::endl;
    std::cout << "M(0,1,1,1,3,0) = " << M(0,1,1,1,3,0) << std::endl;
    std::cout << "M(0,2,1,2,3,1) = " << M(0,2,1,2,3,1) << std::endl;
    std::cout << "M(0,2,1,2,3,0) = " << M(0,2,1,2,3,0) << std::endl;



    return 0;
}