#ifndef NUSSINOV_HPP
#define NUSSINOV_HPP


#include <vector>
#include <string>
#include <utility>

#include "utilities.hpp"
#include "global_variables.hpp"

// Déclaration de la fonction FillMatrix
void FillMatrix(const std::string& sequence, std::vector<std::vector<float>>& m, int theta = ::theta, float pair_energy = ::pair_energy);

float Nussinov(const std::string& sequence,int start_index, int end_index, int theta = ::theta, float pair_energy = ::pair_energy);

// Déclaration de la fonction Backtrack
void Backtrack(int i, int j, const std::vector<std::vector<float>>& m, const std::string& sequence, std::vector<std::pair<int, int>>& S, int theta = ::theta);

#endif // NUSSINOV_HPP