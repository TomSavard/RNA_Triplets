#include "nussinov.hpp"
#include "utilities.hpp"
#include "global_variables.hpp"

// #include <ViennaRNA/PS_rna_plot.h>
#include <iostream>
#include <vector>
#include <string>


void FillMatrix(const std::string& sequence, std::vector<std::vector<int>>& m, int theta = 3, int pair_energy = ::pair_energy) {
    /**
     * @brief Fills the minimum energy matrix 
     * 
     * Detailled description : TO BE DONE
     * 
     * @param sequence : the RNA sequence
     * @param m : the matrix to fill
     * @param theta : the minimum distance between paired bases
     * @param pair_energy : the energy of a pair
     * 
     * @return void
     */
    std::cout << "  Start of FillMatrix" << std::endl;
    int len_seq = sequence.size();

    // Initialize the matrix with zeros.
    m.assign(len_seq, std::vector<int>(len_seq, 0));

    // Fill the matrix
    for (int i = len_seq - 1; i >= 0; i--) {
        for (int j = i + theta + 1; j < len_seq; j++) {
            // Case A : pos i without partner
            m[i][j] = m[i + 1][j];
            // Case B : pos i and j form a base pair
            if (can_pair(sequence[i], sequence[j])) {
                m[i][j] = std::min(m[i][j], m[i + 1][j - 1] + pair_energy);
            }
            // Case C : pairing with another base at index k
            for (int k = i + theta + 1; k < j; k++) {
                if (can_pair(sequence[i], sequence[k])) {
                    m[i][j] = std::min(m[i][j], m[i + 1][k - 1] + m[k + 1][j] + pair_energy);
                }
            }
        }
    }

    std::cout << "  End of FillMatrix" << std::endl;
}


void Backtrack(int i, int j, const std::vector<std::vector<int>>& m, const std::string& sequence, std::vector<std::pair<int, int>>& S, int theta = 3){
    /**
     * @brief Backtrack the minimum energy matrix to find the secondary structure
     * 
     * Detailled description : TO BE DONE
     * 
     * @param i : the first index of the subsequence
     * @param j : the last index of the subsequence
     * @param m : the minimum energy matrix
     * @param sequence : the RNA sequence
     * @param S : the secondary structure
     * @param theta : the minimum distance between paired bases
     * 
     * @return modifies S by reference
     */
    if (j-i < theta) {
        return;
    }
    else {
        // Case A : pos i without partner
        if (m[i][j] == m[i + 1][j]) {
            Backtrack(i + 1, j, m, sequence, S, theta);
            return;
        }
        else if (m[i][j] == m[i+1][j-1] + pair_energy && can_pair(sequence[i], sequence[j])){
            S.push_back(std::make_pair(i,j));
            Backtrack(i + 1, j - 1, m, sequence, S, theta);
            return;
        }
        else {
            for (int k = i + theta + 1; k < j; k++) {
                if (m[i][j] == m[i+1][k-1] + m[k+1][j] + pair_energy && can_pair(sequence[i], sequence[k])) {
                    S.push_back(std::make_pair(i,k));
                    Backtrack(i+1, k - 1, m, sequence, S, theta);
                    Backtrack(k + 1, j, m, sequence, S, theta);
                    break;
                }
            }
        }

    }

}




#ifdef NUSSINOV_TEST
int main() {
    std::string a_sequence = "GCAACUGGCACAAAGGCCUCCUGG";
    int len_seq = a_sequence.size();

    std::cout << "  Input Sequence: " << a_sequence << std::endl;
    std::cout << "  Size of the sequence : " << len_seq << std::endl;

    std::vector<std::vector<int>> m;
    FillMatrix(a_sequence, m);
    print_matrix(m,a_sequence);

    std::vector<std::pair<int, int>> S;
    Backtrack(0, len_seq - 1, m, a_sequence, S, 3);

    std::string dot_bracket_structure = displaySS(S, len_seq);
    std::cout << "  Secondary structure: " << dot_bracket_structure << std::endl; 


    auto [base_pairs, length] = parseSS(dot_bracket_structure);
    std::cout << "  Parsed Base-Pairs: ";
    for (const auto& pair : base_pairs) {
        std::cout << "(" << pair.first << ", " << pair.second << ") ";
    }
    std::cout << std::endl;
    
    // PS_rna_plot(sequence.c_str(), dot_bracket_structure.c_str(), "RNA_plot.ps");

    return 0;
}
#endif //NUSSINOV_TEST