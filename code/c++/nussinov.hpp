#ifndef NUSSINOV_HPP
#define NUSSINOV_HPP


#include <vector>
#include <string>
#include <utility>

#include "utilities.hpp"

// Déclaration de la fonction FillMatrix
void FillMatrix(const std::string& sequence, std::vector<std::vector<int>>& m, int theta = 3, int pair_energy = ::pair_energy);

float Nussinov(const std::string& sequence,int start_index, int end_index, int theta = 3, int pair_energy = ::pair_energy);


// Déclaration de la fonction Backtrack
void Backtrack(int i, int j, const std::vector<std::vector<int>>& m, const std::string& sequence, std::vector<std::pair<int, int>>& S, int theta = 3);

#endif // NUSSINOV_HPP