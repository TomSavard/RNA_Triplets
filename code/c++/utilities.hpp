#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include "global_variables.hpp"

bool can_pair(char base1, char base2);
void print_pairing(const std::string& sequence, const std::vector<std::pair<int,int> >& S);
std::string displaySS(const std::vector<std::pair<int,int> >& S,int size);
std::pair<std::vector<std::pair<int,int> >,int> parseSS(const std::string& structure);
void print_matrix(const std::vector<std::vector<float>>& m, const std::string& sequence, int cellWidth = 3);



#endif //UTILITIES_HPP