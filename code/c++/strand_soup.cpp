#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <ctime>
#include <limits>

#include "nussinov.hpp"
#include "utilities.hpp"

// just for the colors in the terminal
#define RESET   "\033[0m"
#define RED     "\033[31m"
#define GREEN   "\033[32m"



//======================================    The Classes    ==============================================//



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





class output_backtrack{
    /**
     * @brief Class to handle the output of the backtracking
     * 
     * Detailled description : TO BE DONE
     * 
     * @param list_of_sequences : the list of sequences
     * @param list_of_pairs : the list of pairs (x,i,y,j) with x the index of the first strand, i the index of the first nucleotide within x, y the index of the second strand and j the index of the second nucleotide within y.
     * 
     * 
     */
    public:
    void add_sequence(int sequence){
        list_of_sequences.push_back(sequence);
    }
    void add_sequence_front(int sequence){
        list_of_sequences.insert(list_of_sequences.begin(),sequence);
    }
    void add_pair(int x, int i, int y, int j){
        list_of_pairs.push_back({x,i,y,j});
    }
    void shift(int shift){
        for (int i=0; i<int(list_of_pairs.size()); i++){
            list_of_pairs[i][0] += shift;
            list_of_pairs[i][2] += shift;
        }
    }

    void merge(const output_backtrack& sub_output, int shift){// we merge the output of the sub_output with the current output. The merge is done on the right. [upper_output, sub_output_shifted]
        for (const auto& seq : sub_output.list_of_sequences){
            list_of_sequences.push_back(seq + shift);
        }
        for (const auto& pair : sub_output.list_of_pairs){
            list_of_pairs.push_back({pair[0] + shift, pair[1], pair[2] + shift, pair[3]});
        }}

    void print() const{
        std::cout << "Ordered list of sequences: ";
        for (const auto& seq : list_of_sequences){
            std::cout << seq << " ";
        }
        std::cout << std::endl;

        std::cout << "Pairs:  ";
        for (const auto& pair : list_of_pairs){
            std::cout << "(" << pair[0] << "," << pair[1] << "," << pair[2] << "," << pair[3] << ") ";
        }
        std::cout << std::endl;
    }

    private:
    std::vector<int> list_of_sequences;
    std::vector<std::vector<int>> list_of_pairs;
};









//======================================    Energy matrix functions    ==============================================//





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
                // std::cout << "CASE2 with k = " << k << std::endl;
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
                            // std::cout << "CASE3 with k = " << k << std::endl;
                            min_value = value;
                        }
                    }
                }
            }
        }   
    }

    // Case 4 : i is paired to k of strand r

    for (int k=1; k<=j; k++){
        if (can_pair(strands.at(s)[i],strands.at(r)[k])){
            if (k == int(strands.at(r).length()-1)){
                float value = pair_energy + M(m,s,i+1,r,k-1,0);
                if (value < min_value){
                    // std::cout << "CASE4 with k = " << k << std::endl;
                    min_value = value;
                }}
            else {
                float value = pair_energy + M(m,s,i+1,r,k-1,0) + Nussinov_matrices[strands.at(r)](k+1,j);
                if (value < min_value){min_value = value;
                    // std::cout << "CASE4 with k = " << k << std::endl;
                }}
        }
    }
    // if (min_value == M(m,s,i+1,r,j,c)){std::cout << "Case 1" << std::endl;}
    return min_value;
}






void MainAuxiliaryMatrix(std::unordered_map<int, std::string> strands, Matrix6D& M, std::unordered_map<std::string, Matrix2D> Nussinov_matrices){
    /**
     * @brief Computes the 6D matrix of the strand soup.
     * 
     * Detailled description : Sebastian Will suggested the idea of only computing the "square" 6D matrix and replacing the first one by a function with a disjunction of cases.
     * 
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
                                    // std::cout << "CASE GENERAL" << std::endl;
                                    M(m,s,i,r,j,c) = GeneralCaseMinimization(m,s,i,r,j,c,strands,M,Nussinov_matrices);
                                }
                                
                            }
                            //just displaying the values woth colors (red for borders, green for acceptable values (c==1), white for the rest)
                            if (i == 0 || j== 0 || i == M.get_i_size() + 1 || j == M.get_j_size() + 1){
                                std::cout << RED << "M(" << m << "," << s << "," << i << "," << r << "," << j << "," << c << ") = " << M(m,s,i,r,j,c) << RESET << std::endl;
                            }
                            else{
                                if (c ==1){std::cout << GREEN << "M(" << m << "," << s << "," << i << "," << r << "," << j << "," << c << ") = " << M(m,s,i,r,j,c) << RESET << std::endl;}
                                else{std::cout << "M(" << m << "," << s << "," << i << "," << r << "," << j << "," << c << ") = " << M(m,s,i,r,j,c) << std::endl;}
                            }
                        }
                    }
                }
            }
        }
    }
}


 

//======================================    Backtrack functions    ==============================================//




std::vector<int> Find_start_backtrack(Matrix6D& M){
    /**
     * @brief finds the starting point ==> min of the filled energy matrix (apart from the border effect)
     * 
     * Detailled description : TO BE DONE
     * 
     * @param M : the 6D matrix
     * 
     * @return std::vector<int> : the starting point (m,s,i,r,j,c)
     */
    std::vector<int> starting_point(6,-1);
    float min_value = inf_energy;

    for (int s=1; s<=M.get_s_size(); s++){
        for (int r=1; r<= M.get_r_size(); r++){
            if (M(M.get_m_size()-1,s,1,r,M.get_j_size(),1) < min_value){
                min_value = M(M.get_m_size()-1,s,1,r,M.get_j_size(),1);
                starting_point = {M.get_m_size()-1,s,1,r,M.get_j_size(),1};
            }
        }
    }
    if (min_value == inf_energy){
        throw std::runtime_error("No valid starting point found by Find_start_backtrack in the given energy matrix");
    }
    return starting_point;
}





output_backtrack square_backtrack(int m, int s, int i, int r, int j, int c, std::unordered_map<int, std::string> strands, Matrix6D& M, std::unordered_map<std::string, Matrix2D> Nussinov_matrices);





output_backtrack nussinov_backtrack(int s, int i, int j,std::unordered_map<int, std::string> strands, std::unordered_map<std::string, Matrix2D> Nussinov_matrices){
    /**
     * @brief Backtrack the minimum energy matrix to find the secondary structure using the Nussinov algorithm.
     * 
     * Detailled description : TO BE DONE //! we could provide the specific strand and energy matrix instead of everything. It could save space.
     * 
     * @param s : the index of the strand
     * @param i : the first index of the subsequence 
     * @param j : the last index of the subsequence
     * @param strands : the dictionary of strands (int = index, string = sequence)
     * @param Nussinov_matrices : the dictionary of Nussinov matrices (string = sequence, Matrix2D = energy matrix)
     * @param theta : the minimum distance between paired bases 
     * 
     * @return output_backtrack : the secondary structure
     */
    // std::cout << "Nussinov backtrack" << std::endl;
    if (j-i < theta){
        return output_backtrack();
    }
    else {
        // Case A : pos i without partner
        if (Nussinov_matrices[strands.at(s)](i,j) == Nussinov_matrices[strands.at(s)](i+1,j)){
            return nussinov_backtrack(s,i+1,j,strands,Nussinov_matrices);
        }
        // Case B : pos i paired with j
        else if (Nussinov_matrices[strands.at(s)](i,j) == Nussinov_matrices[strands.at(s)](i+1,j-1) + pair_energy && can_pair(strands.at(s)[i],strands.at(s)[j])){
            output_backtrack output = nussinov_backtrack(s,i+1,j-1,strands,Nussinov_matrices);
            output.add_pair(1,i,1,j);
            return output;
        }
        // Case C : pairing with another base at index k
        else {
            for (int k=i+theta+1; k<j; k++){
                if (Nussinov_matrices[strands.at(s)](i,j) == Nussinov_matrices[strands.at(s)](i+1,k-1) + Nussinov_matrices[strands.at(s)](k+1,j) + pair_energy && can_pair(strands.at(s)[i],strands.at(s)[k])){
                    output_backtrack output1 = nussinov_backtrack(s,i+1,k-1,strands,Nussinov_matrices);
                    output_backtrack output2 = nussinov_backtrack(s,k+1,j,strands,Nussinov_matrices);
                    output1.merge(output2,0);
                    output1.add_pair(1,i,1,k);
                    return output1;
                }
            }
        }
    }
    std::cout << "ERROR : Unexpected behaviour in nussinov_backtrack" << std::endl;
    return output_backtrack();
}






output_backtrack bubble_backtrack(int m, int s, int i, int r, int j, int c, std::unordered_map<int, std::string> strands, Matrix6D& M, std::unordered_map<std::string, Matrix2D> Nussinov_matrices){
    /**
     * @brief Backtrack the bubble case.
     * 
     * Detailled description : TO BE DONE
     * 
     * @param m : the number of remaining strands in the soup
     * @param s : the index of the starting strand (1-based)
     * @param i : the index of the starting nucleotide of the first strand (1-based)
     * @param r : the index of the ending strand (1-based)
     * @param j : the index of the ending nucleotide of the last second strand (1-based)
     * @param c : the connectivity bit
     * @param strands : the dictionary of strands (int = index, string = sequence)
     * @param M : the energy_matrix
     * @param Nussinov_matrices : the dictionary of Nussinov matrices (string = sequence, Matrix2D = energy matrix)
     * 
     * @return output_backtrack : the secondary structure
     */

    // We have to backtrack the 6D matrix to find the optimal secondary structure
    // We go from the top to the bottom of the matrix according to the lexico order (m,s,i,r,j,c)

    //CASE 1 : i is left unpaired
    if (M(m,s,i,r,j,c) == M(m,s,i+1,r,j,c)){
        // std::cout << "CASE1" << std::endl;
        return bubble_backtrack(m,s,i+1,r,j,c,strands,M,Nussinov_matrices);
    }

    //CASE 2: i is paired to k of strand s
    for (int k=i+1; k <=int(strands.at(s).length())-1;k++){
        if (can_pair(strands.at(s)[i],strands.at(s)[k])){
            if (M(m,s,i,r,j,c) == pair_energy + Nussinov_matrices[strands.at(s)](i+1,k-1) + M(m,s,k+1,r,j,c)){
                // std::cout << "CASE2" << std::endl;
                output_backtrack output1 = nussinov_backtrack(s,i+1,k-1,strands,Nussinov_matrices);
                output1.add_pair(1,i,1,k); // we add the pair between i and k of strand s 
                output_backtrack output2 = square_backtrack(m,s,k+1,r,j,c,strands,M,Nussinov_matrices);
                output1.merge(output2,0);
                return output1;
            }
        }
    }

    //CASE 3: i is paired to k of new strand t
    if (m>=1){
        for (int t=1; t<=M.get_s_size(); t++){
            for (int m1=0; m1<m; m1++){
                int m2 = m-m1-1;
                for (int k=1; k<=int(strands.at(t).length())-1;k++){
                    if (can_pair(strands.at(s)[i],strands.at(t)[k])){
                        if (M(m,s,i,r,j,c) == pair_energy + M(m1,s,i+1,t,k-1,0) + M(m2,t,k+1,r,j,c)){
                            output_backtrack output1 = square_backtrack(m1,s,i+1,t,k-1,0,strands,M,Nussinov_matrices);
                            output1.add_pair(1,i,1+m1+1,k);
                            output1.add_sequence(t);
                            output_backtrack output2 = square_backtrack(m2,t,k+1,r,j,c,strands,M,Nussinov_matrices); // there is a problem here
                            output2.shift(1+m1);
                            output1.merge(output2,0);
                            return output1;
                        }
                    }
                }
            }
        }
    }

    //CASE 4: i is paired to k of strand r
    for (int k=1; k<=j;k++){
        if (can_pair(strands.at(s)[i],strands.at(r)[k])){
            if (k == int(strands.at(r).length()-1)){
                if (M(m,s,i,r,j,c) == pair_energy + M(m,s,i+1,r,k-1,0)){
                    // std::cout << "CASE4.1" << std::endl;
                    output_backtrack output = square_backtrack(m,s,i+1,r,k-1,0,strands,M,Nussinov_matrices);
                    output.add_pair(1,i,1+m+1,k);
                    return output;
                }
            }
            else{
                if (M(m,s,i,r,j,c) == pair_energy + M(m,s,i+1,r,k-1,0) + Nussinov_matrices[strands.at(r)](k+1,j)){
                    // std::cout << "CASE4.2" << std::endl;
                    output_backtrack output1 = square_backtrack(m,s,i+1,r,k-1,0,strands,M,Nussinov_matrices);
                    output1.add_pair(1,i,1+m+1,k);
                    output_backtrack output2 = nussinov_backtrack(r,k+1,j,strands,Nussinov_matrices);
                    output1.merge(output2,0);
                    return output1;
                }
            }
        }
    }
    std::cout << "ERROR : Unexpected behaviour in bubble_backtrack" << std::endl;
    return output_backtrack();
}


output_backtrack square_backtrack(int m, int s, int i, int r, int j, int c, std::unordered_map<int, std::string> strands, Matrix6D& M, std::unordered_map<std::string, Matrix2D> Nussinov_matrices){
    /**
     * @brief Backtrack the square case.
     * 
     * Detailled description : TO BE DONE
     * 
     * @param m : the number of remaining strands in the soup
     * @param s : the index of the starting strand (1-based)
     * @param i : the index of the starting nucleotide of the first strand (1-based)
     * @param r : the index of the ending strand (1-based)
     * @param j : the index of the ending nucleotide of the last second strand (1-based)
     * @param c : the connectivity bit
     * @param strands : the dictionary of strands (int = index, string = sequence)
     * @param M : the energy_matrix
     * @param Nussinov_matrices : the dictionary of Nussinov matrices (string = sequence, Matrix2D = energy matrix)
     * 
     * @return return 
     */

    // We have to backtrack the 6D matrix to find the optimal secondary structure
    // We go from the top to the bottom of the matrix according to the lexico order (m,s,i,r,j,c)


    if (i>int(strands.at(s).length()-1)){
        if (c==1){
            return output_backtrack(); //nothing
        }
        else {
            if (m==0){
                return nussinov_backtrack(r,1,j,strands,Nussinov_matrices);
            }
            else {
                for (int t=1; t<=M.get_s_size();t++){
                    if (M(m,s,i,r,j,c) == M(m-1,t,1,r,j,1)){
                        output_backtrack output = square_backtrack(m-1,t,1,r,j,1,strands,M,Nussinov_matrices);
                        output.add_sequence(t);
                        output.shift(1);
                        return output;
                    }
                }
            }
        }
    }
    else {
        if (j<1){
            if (c==1){
                return output_backtrack(); //nothing
            }
            else {
                if (m==0){
                    return nussinov_backtrack(s,i,strands.at(s).length()-1,strands,Nussinov_matrices);
                    }
                else {
                    for (int t=1; t<=M.get_s_size();t++){
                        if (M(m,s,i,r,j,c) == M(m-1,s,i,t,strands.at(t).length()-1,1)){
                            output_backtrack output = square_backtrack(m-1,s,i,t,strands.at(t).length()-1,1,strands,M,Nussinov_matrices);
                            output.add_sequence(t);
                            return output;
                        }
                    }
                }
            }
        }
        else{
            return bubble_backtrack(m,s,i,r,j,c,strands,M,Nussinov_matrices);
        }
    }
    std::cout << "ERROR : Unexpected behaviour in square_backtrack" << std::endl;
    return output_backtrack();
}



void full_backtrack( std::unordered_map<int, std::string> strands, Matrix6D& M, std::unordered_map<std::string, Matrix2D> Nussinov_matrices){
    /**
     * @brief Backtrack the 6D matrix to find the optimal secondary structure.
     * 
     * Detailled description : TO BE DONE
     * 
     * @param strands : the dictionary of strands (int = index, string = sequence)
     * @param M : the energy_matrix
     * @param Nussinov_matrices : the dictionary of Nussinov matrices (string = sequence, Matrix2D = energy matrix)
     * 
     * @return void
     */

    std::vector<int> starting_point;
    try{
        starting_point = Find_start_backtrack(M);
    }
    catch (const std::runtime_error& e){
        std::cerr << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }
    std::cout << "Starting point : M(" << starting_point[0] << "," << starting_point[1] << "," << starting_point[2] << "," << starting_point[3] << "," << starting_point[4] << "," << starting_point[5] << ") = " << M(starting_point[0], starting_point[1], starting_point[2], starting_point[3], starting_point[4], starting_point[5]) << std::endl;
    output_backtrack secondary_structure = square_backtrack(starting_point[0],starting_point[1],starting_point[2],starting_point[3],starting_point[4],starting_point[5],strands,M,Nussinov_matrices);
    secondary_structure.add_sequence(starting_point[3]);
    secondary_structure.add_sequence_front(starting_point[1]);
    secondary_structure.print();
}







//=============================================================================================//






int main() {
    std::cout << "\n==================== Strand Soup ====================" << std::endl;
    srand(time(0)); // Initialiser le générateur de nombres aléatoires

    //============= Setting the parameters =============//
    std::cout << "\n========== Setting the parameters ==========" << std::endl;
    int m_start = 3; // Number of sequences to generate
    int sequence_length = 7; // length of the sequences
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
    strands[1] = "$AAAAAAC";
    strands[2] = "$GGGGGGU";
    strands[3] = "$AUUUUUU";
    // strands[3] = "$UUU";
    // strands[4] = "$UUU";
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


    // std::cout << " ========== Test of the class output_backtrack ==========" << std::endl;
    // output_backtrack output;
    // output.print();
    // output.add_sequence(1);
    // output.add_sequence(2);
    // output.print();
    // output.add_pair(1,1,2,2);
    // output.add_pair(1,2,2,3);
    // output.print();
    // output.shift(1);
    // output.print();


    // std::cout << "========== Test of the nussinov_backtrack function ==========" << std::endl;
    // int test_sequence_length = 7;
    // std::unordered_map<int, std::string> test_strands;
    // test_strands[1] = "$AAACUUU";
    // test_strands[2] = "$AAACUUU";
    // std::unordered_map<std::string, Matrix2D> test_Nussinov_matrices;
    // for (const auto& pair : test_strands) {
    //     const auto& seq = pair.second;
    //     std::cout << "  Energy matrix for sequence : " << seq << std::endl;
    //     Matrix2D test_energy_matrix(test_sequence_length, test_sequence_length);
    //     FillMatrix(seq, test_energy_matrix);
    //     test_Nussinov_matrices[seq] = test_energy_matrix;
    //     print_matrix(test_energy_matrix,seq);
    //     std::cout << std::endl;
    // }
    // output_backtrack test_output = nussinov_backtrack(1,1,test_sequence_length,test_strands,test_Nussinov_matrices);
    // test_output.print();


    std::cout << "========== Test of the full_backtrack function ==========" << std::endl;

    std::cout << " ====== Recall ====== " << std::endl;
    std::cout << "  Number of strands : m = " << m_start << std::endl;
    std::cout << "  Length of the sequences : " << sequence_length << std::endl;
    std::cout << "  Theta : " << theta << std::endl;
    std::cout << "  Pair energy : " << pair_energy << std::endl;
    std::cout << "Generated strands:" << std::endl;
    for (const auto& pair : strands){
        std::cout << "  Index " << pair.first << " : " << pair.second << std::endl;
    }
    std::cout << " ==================== " << std::endl;


    full_backtrack(strands, M, Nussinov_matrices);
    // std::vector<int> starting_point = Find_start_backtrack(M);
    // std::cout << "Starting point : M(" << starting_point[0] << "," << starting_point[1] << "," << starting_point[2] << "," << starting_point[3] << "," << starting_point[4] << "," << starting_point[5] << ") = " << M(starting_point[0], starting_point[1], starting_point[2], starting_point[3], starting_point[4], starting_point[5]) << std::endl;

    // output_backtrack secondary_structure = square_backtrack(2,1,1,3,3,1,strands,M,Nussinov_matrices);
    // secondary_structure.add_sequence(starting_point[3]);
    // secondary_structure.add_sequence_front(starting_point[1]);
    // secondary_structure.print();
    return 0;
}