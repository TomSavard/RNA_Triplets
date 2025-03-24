#ifndef NUSSINOV_HPP
#define NUSSINOV_HPP


#include <vector>
#include <string>
#include <utility>

#include "utilities.hpp"
#include "global_variables.hpp"

// Déclaration de la fonction FillMatrix
void FillMatrix(const std::string& sequence, Matrix2D& m, int theta= ::theta, float pair_energy = ::pair_energy);


// Déclaration de la fonction Backtrack
void Backtrack(int i, int j, const Matrix2D& m, const std::string& sequence, std::vector<std::pair<int, int>>& S, int theta = ::theta);

#endif // NUSSINOV_HPP