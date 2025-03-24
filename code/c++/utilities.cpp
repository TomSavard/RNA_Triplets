#include "utilities.hpp"

#include <iostream>
#include <string>
#include <utility>
#include <iomanip>
#include <stack>

Matrix2D::Matrix2D() : rows(0), cols(0) {}
Matrix2D::Matrix2D(int rows, int cols) : rows(rows), cols(cols) {
    data = std::vector<std::vector<float>>(rows, std::vector<float>(cols, 0));}
Matrix2D::Matrix2D(const std::vector<std::vector<float>>& matrix) {
    rows = matrix.size();
    cols = matrix.size();
    data = matrix;
}
float& Matrix2D::operator()(int i, int j) {
    return data[i - 1][j - 1];}
const float& Matrix2D::operator()(int i, int j) const {
    return data[i - 1][j - 1];}
int Matrix2D::get_rows() const {
    return rows;
}
int Matrix2D::get_cols() const {
    return cols;
}


Matrix2D create_matrix2d(const std::vector<std::vector<float>>& data) {
    return Matrix2D(data);
}



bool can_pair(char base1, char base2){ // seems ok
    for (const auto& pair : possible_pairs){
        if (pair.first == base1 && pair.second == base2){
            return true;
        }
    }
    return false;
};



// rappel : l'utilisation de & permet un passage d'argument par référence. On accède à l'objet plutôt que de travailler sur une copie locale. C'est plus économe. On utilise const pour restreindre les modifs.
void print_pairing(const std::string& sequence, const std::vector<std::pair<int,int> >& S){ // seems ok
    for (const auto& pair : S){
        int i = pair.first;
        int j = pair.second;
        if (i>=0 && i<int(sequence.size()) && j>=0 && j<int(sequence.size())){
            std::cout << sequence[i] << " " << sequence[j] << std::endl;
        }
        else {
            std::cerr << "print_pairing - Invalid indices: (" << i << ", " << j << ")" << std::endl;

        }
    }
};





std::string displaySS(const std::vector<std::pair<int,int> >& S,int size){
    std::string structure(size,'.');

    for (const auto& pair : S){
        int i = pair.first;
        int j = pair.second;

        if (i>= 1 && i <= size && j >= 1 && j <= size){
            structure[i-1]='(';
            structure[j-1]=')';
        }
        else {
            std::cerr << "Displayss - Invalid indices: (" << i << ", " << j << ")" << std::endl;
        }
    }
    return structure;
}




std::pair<std::vector<std::pair<int,int> >,int> parseSS(const std::string& structure){
    std::vector<std::pair<int,int>> base_pairs;
    std::stack<int> open_brackets;
    int n = structure.size();
    for (int i=0; i<n; i++){
        if (structure[i] == '('){
            open_brackets.push(i+1);}
        else if (structure[i] == ')'){
            if (open_brackets.empty()){
                std::cerr << "parseSS - Unmatched closing bracket at position " << i+1 << std::endl;
            }
            else {
                int j = open_brackets.top();
                open_brackets.pop();
                base_pairs.push_back(std::make_pair(j,i+1));
            }
        }
    }
    return {base_pairs, n};
};


void print_matrix(const Matrix2D& m, const std::string& sequence, int cellWidth){
    /**
     * @brief print a matrix of floats with a given cell width
     * 
     * @param m (Matrix2D): the matrix to print
     * @param sequence (std::string): the sequence associated with the matrix
     * @param cellWidth (int): the width of the cells. Parameter for the display (default = 3)
     * 
     * @return void
     * 
     */
    int n = m.get_rows();
    std::cout << std::string(cellWidth, ' '); // Leading spaces for alignment
    for (int i = 1; i <= n; ++i) {
        std::cout << std::setw(cellWidth) << sequence[i] << " "; 
    }
    std::cout << std::endl;
    for (int i = 1; i<=n; i++){
        std::cout << std::string(cellWidth, ' ');
        for (int j = 1; j<=n; j++){
            std::cout << std::setw(cellWidth) << m(i,j) << " ";
        }
        std::cout << std::endl;
    }
};




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





#ifdef UTILITIES_TEST
int main() {
    std::string sequence = "GCAACUGGCAC";
    std::vector<std::pair<int, int>> S = {{0, 10}, {1, 9}, {2, 8}};

    std::cout << "  Input Sequence: " << sequence << std::endl;
    std::cout << "  Base-Pairs: ";
    for (const auto& pair : S) {
        std::cout << "(" << pair.first << ", " << pair.second << ") ";
    }
    std::cout << std::endl;


    int size = 11;
    std::cout << "  Size: " << size << std::endl;

    std::string structure = displaySS(S, size);
    std::cout << "  Secondary structure: " << structure << std::endl; 

    auto [base_pairs, length] = parseSS(structure);
    std::cout << "  Parsed Base-Pairs: ";
    for (const auto& pair : base_pairs) {
        std::cout << "(" << pair.first << ", " << pair.second << ") ";
    }
    std::cout << "\n  Length of the string: " << length << std::endl;


    std::cout << "\n Test of the function print_matrix" << std::endl;


    Matrix2D matrix = create_matrix2d({
        {0,0,0,0,1,1,1,1,2,2,2},
        {0,0,0,0,0,0,1,1,1,1,1},
        {0,0,0,0,0,0,0,0,0,1,1},
        {0,0,0,0,0,0,0,0,0,1,1},
        {0,0,0,0,0,0,0,0,0,1,1},
        {0,0,0,0,0,0,0,0,0,1,1},
        {0,0,0,0,0,0,0,0,0,0,1},
        {0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0},
    });

    print_matrix(matrix,sequence);


    std::cout << "\n Test of the function generate_random_sequence" << std::endl;
    std::cout << "  Random sequence: " << generate_random_sequence(10) << std::endl;


    return 0;
}
#endif //UTILITIES_TEST