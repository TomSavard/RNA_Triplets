#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include "global_variables.hpp"

class Matrix2D{
    /**
     * @brief Class for a 2D matrix of floats
     * 
     * The class is used to store a 2D matrix of floats using a 1-based index instead of the 0-based of c++. It provides methods to access and modify the elements. 
     * 
     * @param rows (int): the number of rows
     * @param cols (int): the number of columns
     * 
     * @return void
     * 
     */
    public:
    Matrix2D();
    Matrix2D(int rows, int cols);
    Matrix2D(const std::vector<std::vector<float>>& matrix);

    float& operator()(int i, int j);
    const float& operator()(int i, int j) const;

    int get_rows() const;
    int get_cols() const;

    private:
    int rows, cols;
    std::vector<std::vector<float>> data;

};

Matrix2D create_matrix2d(const std::vector<std::vector<float>>& data);
bool can_pair(char base1, char base2);
void print_pairing(const std::string& sequence, const std::vector<std::pair<int,int> >& S);
std::string displaySS(const std::vector<std::pair<int,int> >& S,int size);
std::pair<std::vector<std::pair<int,int> >,int> parseSS(const std::string& structure);
void print_matrix(const Matrix2D& m, const std::string& sequence, int cellWidth = 3);
std::string generate_random_sequence(int length);



#endif //UTILITIES_HPP